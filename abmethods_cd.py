import os
import gzip
import numpy as np

import scanpy as sc
import anndata as ad
import pandas as pd
import statsmodels.api as sm
import loompy

import h5py
import copy

from scipy.stats.mstats import gmean
from scipy.linalg import helmert

class Patient(object):        

    def __init__(self, data_path:dict, antibodies:list):
        """ data_path structure: {'in_folder':data_path, 'out_folder':output_path, 
        'variants_file_path':variants_file_path, 'genotypes_path':genotypes_path}"""
        
        # general file organization 
        self.path = data_path
        self.ab_map = {name: i for (i, name) in enumerate(antibodies)}
        self.ab_names = antibodies
        self.n_abs = len(self.ab_map)
        self.ab_map['IgG1ctrl'] = self.ab_map['IgG1'] # naming convention changed over the experiments 
        self.ab_map['CD45-1'] = self.ab_map['CD45'] # hashing names
        self.ab_map['CD45-2'] = self.ab_map['CD45'] # hashing
        self.ab_map['CD45-3'] = self.ab_map['CD45'] # hashing
        self.ab_map['CD45-4'] = self.ab_map['CD45'] # hashing
        
        # patient specific data
        self.experiment_list = None
        self.experiment_dict = None
        self.pat_andat = None
        self.pat_andat_raw = None
        self.pat_cell_barcodes_ab = []
    
    @classmethod
    def merge_experiments(cls, exp_nr_list:list, data_path:dict, antibodies:list, min_abs:int, max_IgG1:int, covariables:list, components:int, offset:bool, use_andat:bool=False):
        """Batch loads a set of DAb-seq experiments. Requires a list of DAB_seq experiment names. 
        E.g. [abseq9, abseq13,...] """
        
        # create patient instance
        pat = cls(data_path, antibodies)
        pat.experiment_list = exp_nr_list
        pat.experiment_dict = {}
        
        # load experiment data for antibodies from raw files unless use_andat is set to True
        # then use generated andat files from previous run
        andatas = []
        andatas_raw = []
        for e in exp_nr_list:
            exp = Experiment(e, pat.path, antibodies)
            pat.experiment_dict[e] = exp
            if use_andat:
                try:
                    exp.load_andat()
                except OSError:
                    print('did not find the premade andat file for experiment {}: Check paths and names'.format(e))
                    print('attempting to read from raw files, now...')
                    exp.cells_only()
                    exp.to_anndata(min_abs, max_IgG1)
            else:
                exp.cells_only()
                exp.to_anndata(min_abs, max_IgG1)
            exp.GLM_regression(covariables, components, offset)
            andatas.append(exp.andat_corr)
            andatas_raw.append(exp.andat_raw)
        
        # create a unified patient andat file by concatenating experimental data
        exp_n = []
        for e_ in pat.experiment_list:
            try:
                exp_n.append(str(int(e_[-2:])))
            except ValueError:
                exp_n.append(str(int(e_[-1:])))
                
        pat.pat_andat =  andatas[0].concatenate(andatas[1:], 
                                                   join='outer', 
                                                   batch_key='experiment', 
                                                   batch_categories=exp_n, 
                                                   index_unique='-')
        
        pat.pat_andat_raw =  andatas_raw[0].concatenate(andatas_raw[1:], 
                                                           join='outer', 
                                                           batch_key='experiment', 
                                                           batch_categories=exp_n, 
                                                           index_unique='-')
        
        pat.pat_cell_barcodes_ab.append(pat.pat_andat.obs_names)
        
        return pat
    
    def load_genotypes(self, suffix:str, from_loom:bool = False):
        # load genotyping data from hdf5 compressed file 
        self.filename = os.path.join(self.path['out_folder'], suffix + '_genotype.loom')
        
        if from_loom:
            try:
                self.genotypes = sc.read_loom(self.filename)
            except ValueError:
                print('loom file not found, check paths or try to load from raw files')
            return
        with h5py.File(self.path['genotypes_path'], 'r') as f:
            
            # import hdf5 layers into arrays  
            cell_barcodes = copy.deepcopy([c.decode('utf8') for c in f['CELL_BARCODES']])
            variants = copy.deepcopy([v.decode('utf8') for v in f['VARIANTS']])
            
            # cells with no abs, no genotype
            no_abs = list(set(cell_barcodes) - set(np.array(self.pat_cell_barcodes_ab[0])))
            
            genotypes = pd.DataFrame(np.transpose(f['GT']), index=cell_barcodes, columns=variants).sort_index()
            genotypes.index.name = 'cell_barcode'
            sample_name = ['abseq' + c.split('-')[-1] for c in genotypes.index]
            genotypes['sample_name'] = sample_name
            genotypes.set_index([genotypes.index, 'sample_name'], inplace=True)
            self.genotypes_noAb = genotypes.loc[no_abs] # may have to create loom file for this too, won't be loaded if starting from loom
            genotypes = genotypes.drop(index=no_abs)
            genotypes[genotypes.isnull()] = 3
            
            #adata = ad.AnnData(np.array(genotypes), dtype='int')
            #adata.obs['cell_barcode'] = genotypes.index
            #adata.var['variant_name'] = genotypes.columns
            #adata.filename = os.path.join(self.path['out_folder'], self.filename)
            loompy.create(self.filename ,np.array(genotypes), {'cell_barcode':np.array(genotypes.index)}, {'sample_name':np.array(genotypes.columns)})
            del genotypes

            quality = pd.DataFrame(np.transpose(f['GQ']), index=cell_barcodes, columns=variants).sort_index()
            quality.index.name = 'cell_barcode'
            quality = quality.drop(index=no_abs)
            with loompy.connect(self.filename) as ds:
                ds.layers['quality'] = np.array(quality).astype(int)
            #adata.layers['quality'] = np.array(total_depth).astype(int)
            del quality
            
            total_depth = pd.DataFrame(np.transpose(f['DP']), index=cell_barcodes, columns=variants).sort_index()
            total_depth.index.name = 'cell_barcode'
            total_depth = total_depth.drop(index=no_abs)
            with loompy.connect(self.filename) as ds:
                ds.layers['total_depth'] = np.array(total_depth).astype(int)
            #adata.layers['total_depth'] = np.array(total_depth).astype(int)
            del total_depth
            
            alt_depth = pd.DataFrame(np.transpose(f['AD']), index=cell_barcodes, columns=variants).sort_index()
            alt_depth.index.name = 'cell_barcode'
            alt_depth = alt_depth.drop(index=no_abs)
            with loompy.connect(self.filename) as ds:
                ds.layers['alt_depth'] = np.array(alt_depth).astype(int)
            #adata.layers['alt_depth'] = np.array(alt_depth).astype(int)
            del alt_depth
            
            #adata.write(os.path.join(self.path['out_folder'], self.filename))
            with loompy.connect(self.filename) as ds:
                self.harmonize_Abs(ds)
            #self.genotypes = adata
            # calculate vaf - nan for division by 0
            #vaf = np.divide(alt_depth, total_depth)
        return
    
    def load_variants(self):
        # load variant annotations tsv file
        variant_info = pd.read_csv(self.path['variants_file_path'], sep='\t', header=0, index_col=0, low_memory=False)
        variant_info.index.name = 'variant' 
        self.variant_info = variant_info
        return
    
    def filter_genotype(self, min_alt_depth:int, min_total_depth:int, min_quality:int):
        with loompy.connect(self.filename) as genotypes:
            genotypes[genotypes.layers['quality'][:,:] < min_quality] = 3
            genotypes[genotypes.layers['total_depth'][:,:] < min_total_depth] = 3
            genotypes[((genotypes[:,:] == 1) | (genotypes[:,:] == 2)) & (genotypes.layers['alt_depth'][:,:] < min_alt_depth)] = 3
        return
        
    def harmonize_Abs(self, loompy_handle):
        
        no_gen = set(np.array(self.pat_cell_barcodes_ab[0])) - set([i[0] for i in np.array(loompy_handle.ra['cell_barcode'])])
        matches = np.sort(list(set(np.array(self.pat_cell_barcodes_ab[0])) - no_gen))
        
        self.andat_noGen = self.pat_andat[list(no_gen)]
        self.andat_raw_noGen = self.pat_andat_raw[list(no_gen)]
        self.pat_andat = self.pat_andat[matches]
        self.pat_andat_raw = self.pat_andat_raw[matches]
        return


class Experiment(Patient):
    """Class to deal with antibodies of a DAb-seq experiemtn. Is a child class
    to Patient."""
    
    def __init__(self, name, data_path, antibodies):
        self.name = name
        self.cell_bc = None
        self.ab_dat_umi = None
        self.ab_dat_raw = None
        self.gen_dat = None
        self.amplicon_names = None
        self.experiment_idx = None
        
        # data structures for single-cell data
        self.andat_corr = None
        self.andat_raw = None
        
        try:
            self.ab_names, self.path
        except AttributeError:
            super().__init__(data_path, antibodies)
        
    def load_data(self, QC=False):
        """loads a DAb_seq experiments. Requires a DAB_seq experiment prefix. 
        E.g. abseq9, abseq13,...]
        if set to QC=True all barcode groups including non cells are loaded"""
        
        if QC:
            self.QC_run()
        else:
            self.cells_only()
        return


    def cells_only(self):
        """read the gziped panel read depth table and antibody counts for the called cell BCs"""
        with gzip.open(os.path.join(self.path['in_folder'], self.name + '.cells.tsv.gz'), 'rt') as f:
            pre_amp = [line.strip().split('\t') for line in f]
        pre_amp = np.array(pre_amp)
        self.amplicon_names = pre_amp[0,1:]
        self.cell_bc = {bc : i for (i,bc) in enumerate(pre_amp[1:,0])}
        self.gen_dat = pre_amp[1:,1:].astype(int)


        self.ab_dat_umi = np.zeros([len(self.cell_bc), len(self.ab_names)])
        self.ab_dat_raw = np.zeros([len(self.cell_bc), len(self.ab_names)])
        self.experiment_idx = np.zeros(len(self.cell_bc))
        with gzip.open(os.path.join(self.path['in_folder'], self.name + '_umi_counts.tsv.gz'), 'rt') as f:
            pre_ab = [line.strip().split('\t') for line in f]
        for l in pre_ab:
            try:
                self.ab_dat_umi[self.cell_bc[l[0]], self.ab_map[l[1]]] += int(l[-1])
                self.ab_dat_raw[self.cell_bc[l[0]]] += int(l[2])
                self.experiment_idx[self.cell_bc[l[0]]] = int(l[0][-1])
            except KeyError:
                continue
        return
    
    
    def QC_run(self):
        """read the gziped panel read depth table and antibody counts for the called cell BCs and
        all other BCs. This is mostly useful for quality control purposes and if e.g. new cell
        calling functions should be implemented"""
        with gzip.open(os.path.join(self.path['in_folder'], self.name + '.all.tsv.gz'), 'rt') as f:
            pre_amp = [line.strip().split('\t') for line in f]
        pre_amp = np.array(pre_amp)
        self.amplicon_names = pre_amp[0,1:]
        self.cell_bc = {bc : i for (i,bc) in enumerate(pre_amp[1:,0])}
        self.gen_dat = pre_amp[1:,1:].astype(int)

        self.called_cells = np.ones(len(self.cell_bc))*-1
        with gzip.open(os.path.join(self.path['in_folder'], self.name + '.cells.tsv.gz'), 'rt') as f:
            pre_cells = [line.strip().split('\t') for line in f]
        for l in pre_cells:
            try:
                self.called_cells[self.cell_bc[l[0]]] = 1
            except KeyError:
                continue
        
        self.ab_dat_umi = np.zeros([len(self.cell_bc), len(self.ab_names)])
        self.ab_dat_raw = np.zeros([len(self.cell_bc), len(self.ab_names)])
        self.experiment_idx = np.zeros(len(self.cell_bc))
        with gzip.open(os.path.join(self.path['in_folder'], self.name + '_umi_counts.tsv.gz'), 'rt') as f:
            pre_ab = [line.strip().split('\t') for line in f]
        for l in pre_ab:
            try:
                self.ab_dat_umi[self.cell_bc[l[0]], self.ab_map[l[1]]] += int(l[-1])
                self.ab_dat_raw[self.cell_bc[l[0]], self.ab_map[l[1]]] += int(l[2])
                self.experiment_idx[self.cell_bc[l[0]]] = int(l[0][-1])
            except KeyError:
                continue
        return
    
    def load_andat(self):
        self.andat_raw = sc.read(os.path.join(self.path['out_folder'], self.name + '_raw.h5'))
        return

    def filter_cells(self, min_abs:int=100, max_IgG1:float=5):
        """retains only cells with sum(Ab) >= min_abs reads and IgG1 <= median(IgG1)*max_IgG1"""
        
        abs_pass = self.ab_dat_umi.sum(axis=1) >= min_abs
        isotype_pass = self.ab_dat_umi[:, self.ab_map['IgG1']] <= np.median(self.ab_dat_umi[:, self.ab_map['IgG1']]) * max_IgG1
        
        retain = np.vstack([abs_pass,isotype_pass]).all(axis=0)
        return retain
    
    def to_anndata(self, min_abs:int=100, max_IgG1:float=5.):
        """Bulds a Scanpy AnnData objects for the experiment. I currently do some filtering and 
        batch correction here, but eventually will disentangle those steps. (When we figured out how
        that should be handled)"""

        # do some very basic and non stringent filtering before creating AnnData objects
        # currently does not work with genotype files
        
        retain = self.filter_cells(min_abs, max_IgG1)
        
        #create the AnnData object
        cells = np.array(list(self.cell_bc.keys()))[retain]
        var = pd.DataFrame(index=self.ab_names)
        obs = pd.DataFrame(index=cells)
        adata = ad.AnnData(self.ab_dat_umi.astype('float64')[retain], obs=obs, var=var, dtype='float64')
        
        #add covariates and batch anaotation to AnnData object
        adata.obs['batch'] = self.experiment_idx[retain]
        adata.obs['IgG1'] = self.ab_dat_umi[retain, -1]
        adata.obs['amplicon'] = self.gen_dat[retain].sum(axis=1)
        adata.obs['ab_raw'] = self.ab_dat_raw[retain].sum(axis=1)
        adata.obs['ab_umi'] = self.ab_dat_umi[retain].sum(axis=1)
                           
        self.andat_raw = adata
        
        adata.write(os.path.join(self.path['out_folder'], self.name) + '_raw.h5')
        return
    
    def GLM_regression(self, covariables:list, components:int=1, offset:bool=False):
        """Perform general linear model regression of each Ab vector with provided covariables.
        currently supported covariables are a subset from: [IgG1','amplicon','ab_raw','ab_umi']
        Covariables are SVD transformed and the first N components are retained (default 1)
        
        currently uses linear regression on log transformed data and a Gaussian error model."""
        
        # add pseudo count, log transform, normalize and column center the covariable matrix before SVD 
        covar = np.log(np.array(self.andat_raw.obs[covariables])+1)
        covar = covar / covar.max(axis=0)
        U,s,Vt= np.linalg.svd(covar - covar.mean(axis=0), full_matrices=False)
        
        # use left eigenvectors for regression, correlation is scale invariant
        if offset:
            factors = sm.add_constant(U[:,:components], prepend=True)
        else:
            factors = U[:,:components]
        
        # if some Ab's have all zero counts or were not mesured, the regression cannot be performed on the raw matrix
        # skip columns with all zero in regression
        self.Ab_filter = self.andat_raw.X.sum(axis=0) != 0
        
        fit_val = np.zeros([len(U), len(self.ab_names)])
        residuals = np.zeros([len(U), len(self.ab_names)])
        self.params = np.zeros([len(self.ab_names),factors.shape[1]])
        for i in range(len(self.ab_names)):
            if self.Ab_filter[i]:
                linear_model_result = sm.GLM(np.log(self.andat_raw.X[:,i]+1), factors, family=sm.families.Gaussian()).fit()
                fit_val[:,i] = linear_model_result.fittedvalues
                residuals[:,i] = linear_model_result.resid_response
                self.params[i] = linear_model_result.params
        
        
        self.andat_corr = self.andat_raw.copy()
        self.andat_corr.X = residuals
        
        return

    def compositional_transform(self, add_pseudocount:bool=False):
        """calculated the three Aitchison geometry transforms for the Ab counts.
        - alr uses the IgG1 counts as universal reference
        - ilr contrasts are based on the SVD of clr
        
        if add_pseudocount is set to true add 1 to prevent zero division, otherwise
        cells with zero counts in the denominator create inf and have to be filtered
        before downstream analysis, e.g.:
            clr_filter = np.isfinite(clr).all(axis=1)
            clr = clr[clr_filter,:]"""
        if add_pseudocount:
            self.clr_data = np.log((1+self.andat_raw.X)/gmean(self.andat_raw.X+1, axis=1).reshape(-1,1))
            self.alr_data = np.log((1+self.andat_raw.X)/(self.andat_raw.X[:,-1]+1).reshape(-1,1))
            
            U,s,Vt = np.linalg.svd(self.clr_data, full_matrices=False)
            self.ilr_data = np.dot(U*s,helmert(len(s)).T)
        else:
            self.clr_data = np.log(self.andat_raw.X/gmean(self.andat_raw.X, axis=1).reshape(-1,1))
            self.alr_data = np.log(self.andat_raw.X/(self.andat_raw.X[:,-1].reshape(-1,1)))
            
            finite_clr = np.isfinite(self.clr_data).all(axis=1)
            U,s,Vt = np.linalg.svd(self.clr_data[finite_clr,:], full_matrices=False)
            self.ilr_data = np.dot(U*s,helmert(len(s)).T)
        return

    
    
    
    
def scran_like(array, pool_size=100):
    """CURRENTLY NOT WORKING
    generated pooled cell counts to stabilize statistics and normalize against
    IgG1 count"""
    size_idx = np.argsort(array.sum(axis=1))
    even_size = size_idx[::2]
    odd_size = size_idx[1::2]


    igg_counts = array[np.concatenate([even_size[::-1],odd_size]),:]
    #igg_counts = igg_counts/igg_counts.sum(axis=0)
    pool_size_factors =  np.zeros(array.shape)
    lin_equations = np.zeros([len(pool_size_factors), len(pool_size_factors)])
    
    for i in range(len(pool_size_factors)):
        
        indices = range(i,i+pool_size)
        indices = np.take(np.arange(len(pool_size_factors)),indices, mode='wrap')
        lin_equations[i,indices] = 1
        
        for j in range(array.shape[1]):
            pool_size_factors[i,j] = np.sum(igg_counts[:,j].take(indices, mode='wrap')) /  np.sum(igg_counts[:,-1].take(indices, mode='wrap'))
    restore_idx = np.argsort(np.concatenate([even_size[::-1],odd_size]))
    return pool_size_factors[restore_idx], lin_equations[restore_idx]

    
    """
    pool_size_factor, lin_equations = scran_like(counts)
    
    from scipy.optimize import nnls
    import time

    #system = np.vstack([lin_equations, np.diag(np.ones(len(igg_counts)))*10**-6])

    x_fit = np.zeros(pool_size_factor.shape)
    for i in range(pool_size_factor.shape[1]):
    start = time.time()
    x_fit[:,i] = nnls(lin_equations, pool_size_factor[:,i])[0]
    end = time.time()
    print(end - start)
    
    
    
    """
    
    
    