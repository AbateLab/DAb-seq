'''

dab-seq: single-cell dna genotyping and antibody sequencing
ben demaree 7.9.2019

the main script for pipeline execution

'''

import os
import subprocess
import argparse
import copy
import sys
import time
from multiprocessing.pool import ThreadPool
from multiprocessing import Process

# import functions from external files
import resources
import cell_calling

# option to enable slack messages
slack_enabled = True
slack_token_file = '/home/bdemaree/.slack_token'

def slack_message(message, enabled=False, slack_token_file=None):
    # for posting a notification to a slack channel

    if enabled:
        from slackclient import SlackClient

        with open(slack_token_file) as f:
            token = f.readline().strip()

        channel = 'server-alerts'
        sc = SlackClient(token)

        sc.api_call('chat.postMessage', channel=channel,
                    text=message, username='pipelines',
                    icon_emoji=':adam:')

def wait(processes):
    # waits for processes to finish
    return [process.communicate() for process in processes]

def file_summary(samples, filename='file_summary.txt', cfg=False, to_file=True):
    # displays sample files and class variables

    input_msg = '''
####################################################################################
# INPUT FILE SUMMARY
####################################################################################
    \n'''
    print input_msg

    if to_file:
        f = open(filename, 'w')
        f.write(input_msg)

    for sample in samples:

        s = vars(sample)

        for item in s:
            print item, ': ', s[item]
        print '\n'

        if to_file:
            for item in s:
                if type(s[item]) is list:
                    f.write(item + ': ' + ', '.join(s[item]) + '\n')

                else:
                    f.write(item + ': ' + str(s[item]) + '\n')
            f.write('\n')

    # option to write config file
    if cfg:

        config_msg = '''
####################################################################################
# CONFIG FILE
####################################################################################
        \n'''
        print config_msg

        if to_file:
            f.write(config_msg)

        with open(cfg, 'r') as cfg_file:
            for line in cfg_file:
                print line.strip()
                if to_file:
                    f.write(line)

    # option to write to file
    if to_file:
        f.close()

def get_fastq_names(path_to_fastq, paired=True):
    # gets fastq filenames in a given directory and runs some simple checks
    # assumes files are compressed with .gz extensions

    R1_files = []
    R2_files = []

    for file in os.listdir(path_to_fastq):

        if file.endswith('.fastq.gz'):
            R = file.split('_')[-2]  # R1 or R2

            if R == 'R1':
                R1_files += [path_to_fastq + file]

            elif R == 'R2' and paired:
                R2_files += [path_to_fastq + file]

            else:
                print 'Unexpected filename structure. Exiting...'
                raise SystemExit

    if len(R1_files) != len(R2_files) and paired:
        print 'Unpaired FASTQ files exist! Check input files.'
        raise SystemExit

    R1_files.sort()
    R2_files.sort()

    return R1_files, R2_files

def generate_samples(R1_files, R2_files):
    # store sample filenames in Sample objects

    samples = []        # list of sample objects

    # for experiments with only DNA panels
    if dna_only:
        for i in range(0, len(R1_files)):
            # assign filenames to sample types
            # note: using alphabetization pattern which may not exist in future

            sample_num = i + 1
            sample_name = sample_basename + '-' + str(sample_num)

            # dna panel

            panel_r1 = R1_files[i]
            panel_r2 = R2_files[i]

            # set file locations and append to sample object

            panel_r1_temp = temp_dir + panel_r1.split('.fastq.gz')[0].split('/')[-1] + '_temp.fastq'
            panel_r2_temp = temp_dir + panel_r2.split('.fastq.gz')[0].split('/')[-1] + '_temp.fastq'

            panel_barcodes = barcode_dir + sample_name + '_barcodes_panel.json'

            # antibodies - set filenames to empty string

            ab_r1 = ''
            ab_r2 = ''

            ab_r1_temp = ''
            ab_r2_temp = ''

            ab_barcodes = ''

            ab_reads = ''

            samples.append(resources.TapestriSample(sample_num,
                                                    panel_r1,
                                                    panel_r2,
                                                    panel_r1_temp,
                                                    panel_r2_temp,
                                                    ab_r1,
                                                    ab_r2,
                                                    ab_r1_temp,
                                                    ab_r2_temp,
                                                    panel_barcodes,
                                                    ab_barcodes,
                                                    ab_reads))

    # for experiments with antibody data
    else:
        for i in range(0, len(R1_files) / 2):
            # assign filenames to sample types
            # note: using alphabetization pattern which may not exist in future

            sample_num = i + 1
            sample_name = sample_basename + '-' + str(sample_num)

            # dna panel

            panel_r1 = R1_files[i + len(R1_files) / 2]
            panel_r2 = R2_files[i + len(R1_files) / 2]

            # set file locations and append to sample object

            panel_r1_temp = temp_dir + panel_r1.split('.fastq.gz')[0].split('/')[-1] + '_temp.fastq'
            panel_r2_temp = temp_dir + panel_r2.split('.fastq.gz')[0].split('/')[-1] + '_temp.fastq'

            panel_barcodes = barcode_dir + sample_name + '_barcodes_panel.json'

            # antibodies

            ab_r1 = R1_files[i]
            ab_r2 = R2_files[i]

            ab_r1_temp = temp_dir + ab_r1.split('.fastq.gz')[0].split('/')[-1] + '_temp.fastq'
            ab_r2_temp = temp_dir + ab_r2.split('.fastq.gz')[0].split('/')[-1] + '_temp.fastq'

            ab_barcodes = barcode_dir + sample_name + '_barcodes_abs.json'

            ab_reads = ab_dir + sample_name + '_ab_reads.tsv'

            samples.append(resources.TapestriSample(sample_num,
                                                    panel_r1,
                                                    panel_r2,
                                                    panel_r1_temp,
                                                    panel_r2_temp,
                                                    ab_r1,
                                                    ab_r2,
                                                    ab_r1_temp,
                                                    ab_r2_temp,
                                                    panel_barcodes,
                                                    ab_barcodes,
                                                    ab_reads))

    return samples

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='''
    
    dab-seq: single-cell dna genotyping and antibody sequencing
    ben demaree 2020
    
    input requirements:
    -config file defining file paths and variables (dabseq.cfg)
    -raw fastq files (panel and optionally antibody tags)
    -cell and ab barcode csvs
    -panel bed files
    
    requires the following programs in path:
    -gatk
    -bowtie2
    -itdseek (flt3-calling only)
    -samtools
    -bedtools
    -bcftools
    -cutadapt
    -bbmap
    -snpeff
    
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('sample_name', type=str, help='sample basename')
    parser.add_argument('cfg_file', type=str, help='config filename')
    parser.add_argument('--chem', default='V1', choices=['V1', 'V2'], help='chemistry version (V1 or V2) (default: V1)')
    parser.add_argument('--dna-only', action='store_true', default=False, help='option to run dna panel pipeline only')
    parser.add_argument('--skip-flt3', action='store_true', default=False, help='option to skip FLT3-ITD calling')

    args = parser.parse_args()  # parse arguments

    sample_basename = args.sample_name
    chem = args.chem
    cfg_f = args.cfg_file
    dna_only = args.dna_only
    skip_flt3 = args.skip_flt3

    print '''
    
                                                                                                                                              
                                                bbbbbbbb                                                                                      
DDDDDDDDDDDDD                  AAA              b::::::b                                                                                      
D::::::::::::DDD              A:::A             b::::::b                                                                                      
D:::::::::::::::DD           A:::::A            b::::::b                                                                                      
DDD:::::DDDDD:::::D         A:::::::A            b:::::b                                                                                      
  D:::::D    D:::::D       A:::::::::A           b:::::bbbbbbbbb                         ssssssssss       eeeeeeeeeeee       qqqqqqqqq   qqqqq
  D:::::D     D:::::D     A:::::A:::::A          b::::::::::::::bb                     ss::::::::::s    ee::::::::::::ee    q:::::::::qqq::::q
  D:::::D     D:::::D    A:::::A A:::::A         b::::::::::::::::b                  ss:::::::::::::s  e::::::eeeee:::::ee q:::::::::::::::::q
  D:::::D     D:::::D   A:::::A   A:::::A        b:::::bbbbb:::::::b --------------- s::::::ssss:::::se::::::e     e:::::eq::::::qqqqq::::::qq
  D:::::D     D:::::D  A:::::A     A:::::A       b:::::b    b::::::b -:::::::::::::-  s:::::s  ssssss e:::::::eeeee::::::eq:::::q     q:::::q 
  D:::::D     D:::::D A:::::AAAAAAAAA:::::A      b:::::b     b:::::b ---------------    s::::::s      e:::::::::::::::::e q:::::q     q:::::q 
  D:::::D     D:::::DA:::::::::::::::::::::A     b:::::b     b:::::b                       s::::::s   e::::::eeeeeeeeeee  q:::::q     q:::::q 
  D:::::D    D:::::DA:::::AAAAAAAAAAAAA:::::A    b:::::b     b:::::b                 ssssss   s:::::s e:::::::e           q::::::q    q:::::q 
DDD:::::DDDDD:::::DA:::::A             A:::::A   b:::::bbbbbb::::::b                 s:::::ssss::::::se::::::::e          q:::::::qqqqq:::::q 
D:::::::::::::::DDA:::::A               A:::::A  b::::::::::::::::b                  s::::::::::::::s  e::::::::eeeeeeee   q::::::::::::::::q 
D::::::::::::DDD A:::::A                 A:::::A b:::::::::::::::b                    s:::::::::::ss    ee:::::::::::::e    qq::::::::::::::q 
DDDDDDDDDDDDD   AAAAAAA                   AAAAAAAbbbbbbbbbbbbbbbb                      sssssssssss        eeeeeeeeeeeeee      qqqqqqqq::::::q 
                                                                                                                                      q:::::q 
                                                                                                                                      q:::::q 
                                                                                                                                     q:::::::q
                                                                                                                                     q:::::::q
                                                                                                                                     q:::::::q
                                                                                                                                     qqqqqqqqq
                                                                                                                                              
                                  `,;*########+*:.                                                                                                   
                               `,*##+i::,,,,,:::;*##+,                                                                                                
                            `:+#+**:,,,,,,,,,:**:,,,;+#;                                                                                              
                          .*#*:,+nni,,,,,,,,,;z#:,,,,,:*#.                                                                                            
                        :#+;,,,,*#*:,,,,,,,,,,::,,,,,,,,:#;                                                                                           
                      ;#*:,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:+*                                                                                          
                    :#*:,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,+*                                                                                         
                  ,zn;,,,,,,,,,,,,,,,,,,:ii:,,,,,,,,,,,,,,,,+*                                                                                        
                `*xnz;,,,,,,,,,,,,,,,,,:#nn+,,,,,,,,;i,,,,,,,+;                                                                                       
               ,#i##;,,,,,:i*i:,,,,,,,,:+#+:,,,,,,,:zz;,,,,,,:z,                                                                                      
              i#:,,,,,,,,:#nnz:,,,,,,,,,,:,,,,,,,,,:+*:,,,,,,,:z`                                                                                     
            `#i,,,,,,,,,,;znzi,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,i*                                                                                     
           .#;,,,,,,,,,,,,;;:,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#.                                                                                    
          ,#:,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,;+                                                                                    
         :#:,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#`                                                                                   
        :#:,,,,::,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,ii                                                                                   
       ,+:,,,,;#+,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:#                                                                                   
      .n;,,,,:#n#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,z`                                                                                  
     `nn;,,,,;z#:,,,,,,,:;*+#*:,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#,                                                                                  
     +n#:,,,,,;:,,,,,,;+z#*;::,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*;                                                                                  
    :+i:,,,,,,,,,,,,;##;,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,;i                                                                              ;;: 
   `#:,,,,,,,,,,,,:+#;,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,;*                                                                           :++;,i,
   ii,,,,,,,,,,,,:z*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:+                                                                     .,;*+#*:,,,:*
  `#,,,,,,,,,,,,;z;,,,,,,,,,,,:;ii,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:+                                                             .:;*++##+*i::,,,,,,,+
  ii,,,,,,,,,,,:z:,,,,,,,,:i#zz#*;,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,;*                                                        .;*##+*i:::,,,,,,,,,,,,,,+
  #,,,,,,,,,,,:z;,,,,,,,;#z+;:,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,;*                                                      i#+i:,,,,,,,,,,,,,,,,,,,,,:*
 .+,,,,**,,,,,+i,,,,,,:##;,,,,,,,,,,,,,,,,:::,,,,,,,,,,,,,,,,,,,,,ii                                                   `*#;:,,,,,,,,,,,,,,,,,,,,,,,:;;
 *;,,,inz:,,,:#,,,,,,*z;,,,,,,,,,,,,,,:i#zz#i,,,,,,,,,,,,,,,,,,,,,+,                                                  .#;,,,,,,,,,,,,,,,,,,,,,,,,:;;#`
 #:,,,#n#,,,,;;,,,,:#+,,,,,,,,,,,,:i#z#*;:,,,,,,,,,,,,,,,,,,,,,,,,z`                                                 ,#:,,,,,,,,,,,,,,,,,,,,,,,,:;;** 
`#,,,:zz;,,,,,,,,,:zi,,,,,,,,,,,;+z+;:,,,,,,,,,,,,,,,,,,,,,,,,,,,:#                                                 ,#:,,,,,,,,,,,,,,,,,,,,,,,,:;;*z` 
:*,,,,i;,,,,,,,,,:z;,,,,,,,,,:iz#;,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*;                                                `#:,,,,,,,,,,,,,,,,,,,,,,,:;;;zM,  
i;,,,,,,,,,,,,,,:z;,,,,,,,,,*z+:,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,z`                                                +;,,,,,,,,,,,,,,,,,,,,,,,:;;;#Mi   
*:,,,,,,,,,,,,,,#i,,,,,,,,iz*:,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,;*                                                ;*,,,,,,,,,,,,,,,,,,,,,,:;;;;*Wi    
#*,,,,,,,,,,,,,i+,,,,,,,:##:,,,,,,,;+####+;,,,,,,,,,,,,,,,,,,,,,#.                                                #:,,,,,,,,,,,,,,,,,,,,:;ii;;;zi     
zz,,,,,,,,,,,,:#:,,,,,,iz;,,,,,,,:#+:,,,,;+#:,,,,,,,,,,,,,,,,,,i*                                                ,*,,,,,,,,,,,,,,,,,,::;+nxx*;#;      
zz,,,,,,,,,,,,*i,,,,,,*#:,,,,,,,,#;........:#*,,,,,,,,,,,,,,,,:z`                                               .z:,,,,,,,,,,,,,,::;;;inn##x*z:       
#*,,,,,,,,,,,,i,,,,,,++,,,,,,,,,i+..........,++:,,,,,,,,,,,,,,+;                                               :#+,,,,,,,,,,,,,::;;;;;zz+#xn#.        
+;,,,,,,,,,,,,,,,,,,+*,,,,,,,,,,#,............i#,,,,,,,,,,,,,;#                                              `**;:,,,,,,,,,,,,:;;;;;;;xnnxn;          
ii,,:i:,,,,,,,,,,,,+*,,,,,,,,,,:z..............i+,,,,,,,,,,,:#.                                             ,#;,,,,,,,,,,,,,,:;;;;;;;;i*z+`           
:*,,ini,,,,,,,,,,,i+,,,,,,,,,,,;#...............**,,,,,,,,,,*i                                            `*+:,,,,,,,,,,,,,,:;;;;*#z#*#+.             
.#,,+n*,,,,,,,,,,:z:,,,,,,,,,,,;#................#;,,,,,,,,;#                                            ,#i,,,,,,,,,,,,,,,:;;;inn#zx#,               
`#,,+ni,,,,,,,,,,#;,,,,,,,,,,,,:#................,z:,,,,,,iz.                                          `*+:,,,,,,,,,,,,,,,,:;;;nzzn+.                 
 #:,;i:,,,,,,,,,;#,,,,,,,::::,,,#.................i*,,,,,,**                                          ,#i,,,,,,,,,,,,,:::,,:;;+Mzi`                   
 ii,,,,,,,,,,,,,#:,,,,,;##++##i:#,......,:.........#:,,,,,:#                                        `*+:,,,,,,,,,,,,,;####++++i.                      
 .#,,,,,,,,,,,,:#,,,,,*#,....,*z#i.....*nx;........;+,,,,,,*;   `:i++++i.                          :#;,,,,,,,,,,,,,;#*.                               
  #:,,,,,,,,,,,i*,,,,;#,.......,+n....,nxxz,.......,#:,,,,,,+.i++i;:,,,;++`                      `++:,,,,,,,,,,,,i#i`                                 
  ;*,,,,,,,,,,,+;,,,,#:..........#;....zxxxi........ii,,,,,:##;,,,,,,,,,,;#`                    ;#;,,,,,,,,,,,,i#i`                                   
  `#,,,,,,,,,,,::,,,:z...........:#....ixxxn,.......,#,,,,iz;,,,,,,,,,,,,,:#`                 .+*,,,,,,,,,,,,i#i`                                     
   ii,,,,,,,,,,,,,,,:+............+:...,nxxx#,.......#:,,:i,,,,,,,,,,::::,,i;                ;+:,,,,,,,,,,:i#;`                                       
   `#:,,,,,,,,,,,,,,:+............,#,...*xxxxi.......+;,,,,,,,,,,:i#####z#ii;              .+i,,,,,,,,,,,i#i`                                         
    :+,,,,,,,,,,,,,,:+.............*i...,zxxxxi......ii,,,,,,,,;##*:,,,,,:*n#,`           i+:,,,,,,,,,,izi`                                           
     +i,,,,,,,,,,,,,:#.............,#,...:nxxxx;.....ii,,,,,,:##;,,,,,,,,,,:;+##*;.`    .#i,,,,,,,,,,i#*`                                             
     `#;,,,,,,,,,,,,,z......:;,.....;+....;xxxxx:....+;,,,,,iz;,,,,,,,,,,,,,,,,::i+##*;i#:,,,,,,,,,;#i`                                               
      `#;,,,,,,,,,,,,#:....*xx*......+i....;nxxx:...,z,,,,:#+:,,,,,,,,,,,,,,,,,,,,,,,:i+z;,,,,,,,;#*`                                                 
       `+*:,,,,,,,,,,i*...,xxxx:.....,#;....:zx*....i*,,,:zi,,,,,,,,,,,,,,,,,,,,,,,,,,,,:z,,,,,;#*`                                                   
         ;#i,,,,,,,,,:z,...#xxxz,....ii#*....,:....:z:,,;z;,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#,,,;#*`                                                     
          `i#+i;:::;iizi...:nxxx+,...#:,*#:.......:z:,,*#:,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:#,;#*`                                                       
             .;*++*i;:,#,...*xxxx*..,#,,,:##i,.,:+#:,,++:,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:#*#*`                                                         
                       ;*...,zxxxx*.:+,,,,,:*#nz+;,,:#*,,,,,,,,,:::;::::::::,,,,,,,,::;zM*`                                                           
                        +:...:nxxxx#ii,,,,,,,,i+,,,:z;,,,,,,,,:;;;;;;;;i**i;;;;;;;;;;in+.                                                             
                        `#,...;nxxxxx;,,,,,,,,,z;,;z;,,,,,,,,;;;;;;;;;zxnxxni;;;i#nz##,                                                               
                         .#,...:nxxxM:,,,,,,,,,:z*z:,,,,,,;+;;;iii;;;*n##++z+;;ixznz,                                                                 
                          ,#,...,#xnx:,,,,,,,,,,*z:,,,,,,:z,z*znnnni;*x###nxi;;zxz;                                                                   
                           .#:....;,#:,,,,,,,,,*+,,,,,,,:#, `*x##+#xi;*zzz+ii+nx,                                                                     
                            `+*,....z:,,,,,,,,++,,,,,,,,+:    ,xMnnM+*+++#znnn#n:                                                                     
                             ,zz#++#x:,,,,,,:#*,,,,,,,,*nn+` ;xz##zznnnzz#######i                                                                     
                            :#:,:;;:z:,,,,,:#i,,,,,,,,in##zn#n##################*                                                                     
                           .#:,,,,,,z:,,,,:zi,,,,,,,,;xx####Mz#################zi                                                                     
                          `#:,,,,,,,+;,,,:z;,,,,,,,,:z:*x####M#################n:                                                                     
                          ;i,,,,,,,,**,,:z;,,,,,,,,:z;,,;n###zx##############n#.                                                                      
                         `#,,,,,,,,;##,:z;,,,,,,,,,#i,,,,:n###x############zn:                                                                        
                         :*,,,,,,:+#;z:z;,,,,,,,,,+*,,,,,,;n##zn#########nn;`                                                                         
                         +:,,,,,:#i,,#z:,,,,,,,,,i+,,,,,,,,z###x#######nz;`                                                                           
                         #,,,,,:z;,,:z;,,,,,,,,,;#:,,,,,,,,;Mn#x###n#nz:                                                                              
                        `#,,,,:z:,,:z;,,,,,,,,,:z:,,,,,,,,,,nzxM###zM+                                                                                
                        `#,,,,#;,,:z;,,,,,,,,,,#;,,,,,,,,,,,z##n#####ni                                                                               
                        `#,,,i*,,:#;,,,,,,,,,,**,,,,,,,,,,,,+#########ni                                                                              
                         #,,,:,,,+i,,,,,,,,,,;#,,,,,,,,,,,,,+##########n:                                                                             
                         *:,,,,,**,,,,,,,,,,:#:,,,,,,,,,,,,,############n.                                                                            
                         ,+,,,,;#,,,,,,,,,,,+i,,,,,,,,,,,,,,n#############                                                                            
                          i#i;iz:,,,,,,,,,,;#,,,,,,,,,,,,,,;n############n:                                                                           
                           .;*ni,,,,,,,,,,:#:,,,,,,,,,,,,,,z##############n`                                                                          
                             .+,,,,,,,,,,,**,,,,,,,,,,,,,,*n##############z*                                                                          
                             +:,,,,,,,,,,:nz:,,,,,,,,,,,,+n################n.                                                                         
                            .+,,,,,,,,,,,##zz;,,,,,,,,,;zn##################+                                                                         
                            *:,,,,,,,,,,;x###nz+;::::i#xz###################n,                                                                        
                           `+,,,,,,,,,,:zz#####znnnxn;.,z####################+                                                                        
                           :i,,,,,,,,,:#z########zni`   #####################x,                                                                       
                           *:,,,,,,,,,+n########n+`     :n####################+                                                                       
                           +:,,,,,,,:+n#######nz,       `n####################n`                                                                      
                           zz:,,,,,inz######nz:          #####################z;                                                                      
                           znn+ii+nn#####zn#,            i######################                                                                      
                           ,zzznnz####znz*.              :z####################n`                                                                     
                            `*znnnnnz+;.                 ,n####################n:                                                                     
                               ....                      .n#####################*                                                                     
                                                         .n#####################z                                                                     
                                                         .n#####################n`                                                                    
                                                         :z#####################n,                                                                    
                                                         *z#####################z;                                                                    
                                                         ########################+                                                                    
                                                        `n#######################z                                                                    
                                                        ,n#######################n                                                                    
                                                        *########################n`                                                                   
                                                        n########################x,                                                                   
                                                       ,n########################n,                                                                   

Initializing pipeline...

    '''

    # load config file variables
    # be careful about using exec
    if not os.path.isfile(cfg_f):
        print 'Config file not found! Please check the file name and path.'
        raise SystemExit

    else:
        with open(cfg_f, 'r') as cfg:
            for line in cfg:
                if line[0] == '#' or line[0] == ' ':
                    continue
                elif line[0] == '[':
                    if 'Antibody' in line and dna_only:
                        line = cfg.next()
                        while line[0] != '[':
                            try:
                                line = cfg.next()
                            except StopIteration:
                                break
                else:
                    var = line.split("#", 1)[0].strip()  # to remove inline comments
                    exec(var)

    # check all files exist
    all_vars = copy.copy(globals())
    input_files = [all_vars[f] for f in all_vars if '_file' in f and f != '__file__']
    input_files.append(bt2_ref + '.1.bt2')
    missing_files = []
    for f in input_files:
        if not os.path.exists(f):
            missing_files.append(f)

    # print missing files, if any, and exit
    if missing_files != []:
        print 'The following input files could not be found:'
        for f in missing_files:
            print f
        print '\nExiting...\n'
        raise SystemExit

    # check that the input fastq directories exist and are not empty
    if not os.path.exists(panel_fastq_dir):
        print 'FASTQ input directory (%s) does not exist! Exiting...\n' % panel_fastq_dir
        raise SystemExit

    elif os.listdir(panel_fastq_dir) == []:
        print 'FASTQ input directory (%s) is empty! Exiting...\n' % panel_fastq_dir
        raise SystemExit

    if not dna_only:
        if not os.path.exists(ab_fastq_dir):
            print 'FASTQ input directory (%s) does not exist! Exiting...\n' % ab_fastq_dir
            raise SystemExit

        elif os.listdir(ab_fastq_dir) == []:
            print 'FASTQ input directory (%s) is empty! Exiting...\n' % ab_fastq_dir
            raise SystemExit

    # create all other directories for this run
    # if it already exists, ignore and continue
    dirs = [all_vars[d] for d in all_vars if '_dir' in d]
    dirs.sort()
    for d in dirs:
        if not os.path.exists(d):
            os.mkdir(d)

    # setup barcode extraction parameters based on chemistry
    if chem == 'V1':
        cell_barcode_csv_file = sys.path[0] + '/barcodes/mb_cell_barcodes_v1.csv'

        r1_start = 'CGATGACG'
        r1_end = 'CTGTCTCTTATACACATCT'
        r2_end = 'CGTCATCG'

        bar_ind_1, bar_ind_2 = range(8), range(-8, 0)

    elif chem == 'V2':
        cell_barcode_csv_file = sys.path[0] + '/barcodes/mb_cell_barcodes_v2.csv'

        r1_start = 'GTACTCGCAGTAGTC'
        r1_end = 'CTGTCTCTTATACACATCT'
        r2_end = 'GACTACTGCGAGTAC'

        bar_ind_1, bar_ind_2 = range(9), range(-9, 0)

    # set minimum read lengths (should rarely need to be modified)
    r1_min_len_panel = 30
    r2_min_len_panel = 25
    r1_min_len_ab = 0
    r2_min_len_ab = 30

    # send slack notification
    start_time = time.time()
    start_time_fmt = str(time.strftime('%m-%d-%Y %H:%M:%S', time.localtime(start_time)))
    slack_message('Pipeline started for sample %s at %s.' % (sample_basename, start_time_fmt),
                  slack_enabled,
                  slack_token_file)

    print '''
####################################################################################
# Step 1: get input file names and store in TapestriSample objects
####################################################################################
    '''

    # get fastq files for dna panel
    R1_files, R2_files = get_fastq_names(panel_fastq_dir)

    # if experiment includes antibody data
    if not dna_only:

        # get fastq files for antibodies
        R1_files_ab, R2_files_ab = get_fastq_names(ab_fastq_dir)

        R1_files += R1_files_ab
        R2_files += R2_files_ab

    R1_files.sort()
    R2_files.sort()

    # store sample info in Sample objects
    samples = generate_samples(R1_files, R2_files)

    # display and write sample summary file
    file_summary(samples, summary, cfg=cfg_f, to_file=True)

    print '''
####################################################################################
# Step 2: filter reads for cell barcode and perform error correction
####################################################################################
'''

    # load mission bio barcode csv file
    barcodes = resources.load_barcodes(cell_barcode_csv_file, 1, False)

    # generate hamming dictionary for error correction
    barcodes = resources.generate_hamming_dict(barcodes)

    print '%d unique barcode sequences loaded into dictionary.\n' % len(barcodes)

    # for panel reads, filter reads with valid barcode structure and export to new fastq
    print 'Extracting barcodes from raw fastq files...\n'

    # first, identify reads with valid barcode

    process_barcodes = []

    for sample in samples:
        # panel files
        p = Process(
            target=sample.filter_valid_reads,
            args=(r1_start,
                  barcodes,
                  bar_ind_1,
                  bar_ind_2,
                  'panel'))
        process_barcodes.append(p)
        p.start()

        if not dna_only:
            # ab files
            p = Process(
                target=sample.filter_valid_reads,
                args=(r1_start,
                      barcodes,
                      bar_ind_1,
                      bar_ind_2,
                      'ab'))
            process_barcodes.append(p)
            p.start()

    # wait for processes to finish
    for p in process_barcodes:
        p.join()

    # second, cut adapters from reads and add barcodes to header

    cut_adapters = []

    for sample in samples:
        # panel files
        p = Process(
            target=sample.barcode_reads,
            args=(r1_start,
                  r1_end,
                  r2_end,
                  r1_min_len_panel,
                  r2_min_len_panel,
                  'panel'))
        cut_adapters.append(p)
        p.start()

        if not dna_only:
            # ab files
            p = Process(
                target=sample.barcode_reads,
                args=(r1_start,
                      r1_end,
                      r2_end,
                      r1_min_len_ab,
                      r2_min_len_ab,
                      'ab'))
            cut_adapters.append(p)
            p.start()

    # wait for processes to finish
    for p in cut_adapters:
        p.join()

    print '''
# ###################################################################################
# Step 3: import ab reads, error correct, and count
# ###################################################################################
'''

    if not dna_only:

        # load ab barcode csv file (with descriptions)
        barcodes = resources.load_barcodes(ab_barcode_csv_file, 1, False)
        barcode_descriptions = copy.deepcopy(barcodes)

        # generate hamming dictionary for error correction
        barcodes = resources.generate_hamming_dict(barcodes)

        # process ab reads and look for barcode
        ab_process = []
        for sample in samples:
            # ab files
            p = Process(
                target=sample.process_abs,
                args=(barcodes,
                      barcode_descriptions,
                      ab_handles,
                      ab_bar_coord,
                      ab_umi_coord,
                      min_umi_qual))
            ab_process.append(p)
            p.start()

        # wait for processes to finish
        for p in ab_process:
            p.join()

        # merge ab reads files into a single file and remove old files
        ab_files = [s.ab_reads for s in samples]

        # concatenate files and remove old ones
        wait([subprocess.Popen('cat %s > %s' % (' '.join(ab_files), ab_reads_merged), shell=True)])
        for f in ab_files:
            try:
                os.remove(f)
            except OSError:
                pass

        # count the ab read umis
        resources.count_umis(ab_reads_merged, umi_counts_merged)

    print '''
###################################################################################
Step 4: count read alignments to inserts
###################################################################################
'''

    # get r1 reads for all panel samples
    panel_r1_files = [s.panel_r1_temp for s in samples]

    # align r1 reads to inserts to obtain read counts across all barcodes
    all_tsv = barcode_dir + sample_basename + '.all.tsv'  # alignment counts for all barcodes
    resources.count_alignments(panel_r1_files, amplicon_file, human_fasta_file, all_tsv, temp_dir)

    print '''
###################################################################################
Step 5: call valid cells using selected method
###################################################################################
'''

    # call valid cells using cell_caller function
    valid_cells = cell_calling.call(barcode_dir, sample_basename, 'second_derivative', False)

    # create SingleCell objects for each valid cell
    cells = [resources.SingleCell(barcode,
                                  by_cell_fastq_dir,
                                  by_cell_bam_dir,
                                  by_cell_gvcf_dir,
                                  by_cell_flt3_dir)
             for barcode in valid_cells]

    print '%s valid cells found!\n' % len(cells)

    print '''
###################################################################################
Step 6: write valid cells from panel reads to separate fastq files
###################################################################################
'''

    # check that all barcodes have the same length
    # this ensures demultiplexing will work correctly
    bar_length = list(set([len(k) for k in valid_cells]))
    assert len(bar_length) == 1, 'All barcodes must have the same length!'

    # split cell fastq files for panel reads associated with valid cells
    barcode_names = temp_dir + 'barcodes.txt'
    with open(barcode_names, 'w') as f:
        for b in valid_cells:
            f.write(b + '\n')

    # split files by cell barcode using bbmap demuxbyname.sh
    split_files = []
    for sample in samples:
        demux_cmd = 'demuxbyname.sh prefixmode=f length=%d in1=%s in2=%s out=%s names=%s' % (
            bar_length[0],
            sample.panel_r1_temp,
            sample.panel_r2_temp,
            by_cell_fastq_dir + '%.fastq',
            barcode_names
        )
        p = subprocess.Popen(demux_cmd, shell=True)
        split_files.append(p)

    # wait for processes to finish
    wait(split_files)

    print '''
####################################################################################
# Step 7: align, convert, sort, index panel reads (optional: call FLT3-ITDs)
####################################################################################
'''

    # limit number of cells to preprocess at a time (based on hardware limitations)
    n_preprocess = 300

    # create pool of workers and run through all samples
    preprocess_pool = ThreadPool(processes=n_preprocess)

    # align and index cells
    for c in cells:
        preprocess_pool.apply_async(resources.SingleCell.align_and_index, args=(c, bt2_ref,))

    preprocess_pool.close()
    preprocess_pool.join()

    # optionally, call FLT3-ITDs using ITDseek
    if not skip_flt3:

        preprocess_pool = ThreadPool(processes=n_preprocess)

        for c in cells:
            preprocess_pool.apply_async(resources.SingleCell.call_flt3, args=(c, human_fasta_file,))

        preprocess_pool.close()
        preprocess_pool.join()

    else:
        os.rmdir(by_cell_flt3_dir)

    print '''
####################################################################################
# Step 8: perform variant calling for all cells
####################################################################################
'''

    # limit number of cells to call variants at a time (based on hardware limitations)
    n_call_variants = 70

    # create pool of workers and run through all samples
    call_variants_pool = ThreadPool(processes=n_call_variants)

    for c in cells:
        call_variants_pool.apply_async(resources.SingleCell.call_variants, args=(c,
                                                                                 human_fasta_file,
                                                                                 interval_file,))

    call_variants_pool.close()
    call_variants_pool.join()

    ################################################################################

    # delete temporary files, if selected
    if fastq_cleanup:
        for s in samples:
            if os.path.isfile(s.panel_r1_temp):
                os.remove(s.panel_r1_temp)
            if os.path.isfile(s.panel_r2_temp):
                os.remove(s.panel_r2_temp)
            if os.path.isfile(s.ab_r1_temp):
                os.remove(s.ab_r1_temp)
            if os.path.isfile(s.ab_r2_temp):
                os.remove(s.ab_r2_temp)

    print 'Pipeline complete!'

    # send slack notification
    elapsed_time = time.time() - start_time
    elapsed_time_fmt = str(time.strftime('%Hh %Mm %Ss', time.gmtime(elapsed_time)))
    slack_message('Pipeline complete for sample %s! Total elapsed time is %s.' % (sample_basename, elapsed_time_fmt),
                  slack_enabled,
                  slack_token_file)