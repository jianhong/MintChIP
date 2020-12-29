#!/usr/bin/env python

#######################################################################
#######################################################################
## Created on Nov. 10, 2020 to create UCSC trackhub file from file list
## Copyright (c) 2020 Jianhong Ou (jianhong.ou@gmail.com)
#######################################################################
#######################################################################

import os
import sys
import glob
import errno
import argparse
import trackhub

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

def parse_args(args=None):
    Description = 'Create UCSC trackhub file from a list of files and associated colours - ".bed", ".narrowPeak", ".broadPeak", ".bw", ".bigwig" files currently supported.'
    Epilog = """Example usage: python create_trackhub.py <OUTPUT_FOLDER> <LIST_FILE> <GENOME>"""

    argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    
    ## REQUIRED PARAMETERS
    argParser.add_argument('OUTPUT_FOLDER', help="Folder for UCSC trackhub files")
    argParser.add_argument('LIST_FILE', help="Tab-delimited file containing two columns i.e. file_name\tcolour. Header isnt required.")
    argParser.add_argument('GENOME', help="Full path to genome fasta file or shorthand for genome available in UCSC e.g. hg19.")
    argParser.add_argument('EMAIL', help="email address")
    argParser.add_argument('DESIGN_FILE', help="design file")
    argParser.add_argument("POSTFIX", help="Postfix of the bigWig file")
    
    ## OPTIONAL PARAMETERS
    argParser.add_argument('-pp', '--path_prefix', type=str, dest="PATH_PREFIX", default='', help="Path prefix to be added at beginning of all files in input list file.")

    return argParser.parse_args(args)


############################################
############################################
## HELPER FUNCTIONS
############################################
############################################

def makedir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

############################################
############################################
## MAIN FUNCTION
############################################
############################################
tracktypes = ['bigWig', 'bam', 'bigBed', 'vcfTabix', 'bigNarrowPeak',
              'bigBarChart', 'bigChain', 'bigGenePred', 'bigBroadPeak',
              'bigMaf', 'bigPsl', 'halSnake']
TrackType = {'bw':'bigWig', 'bb':'bigBed', 'bed':'bigBed', 
             'narrowpeak':'bigNarrowPeak', 'broadpeak':'bigBroadPeak'}
Visibility = {'bw':'full', 'bb':'dense', 'bed':'dense',
             'narrowpeak':'dense', 'broadpeak':'dense'}
for tt in tracktypes:
    TrackType[tt.lower()] = tt
    if tt in ['bam', 'bigBed', 'vcfTabix', 'bigNarrowPeak',
              'bigBarChart', 'bigChain', 'bigGenePred', 'bigBroadPeak',
              'bigMaf', 'bigPsl', 'halSnake']:
      Visibility[tt.lower()] = 'dense'
    else:
      Visibility[tt.lower()] = 'full'

def create_trackhub(OutFolder,ListFile,Genome,EMAIL,DesignFile,Postfix,PathPrefix=''):
    ERROR_STR = 'ERROR: Please check design file'
    HEADER = ['group', 'replicate', 'fastq_1', 'fastq_2']

    makedir(OutFolder)
    
    dIn = open(DesignFile, 'r')
    header = dIn.readline().strip().split(',')
    if header[:4] != HEADER:
        print("{} header: {} != {}".format(ERROR_STR,','.join(header),','.join(HEADER)))
        sys.exit(1)

    paramColn = {}
    for i in range(len(header)):
      if header[i][:6]=="track_": # header start with track_
        paramColn[header[i][6:]]=i

    sampleDesignDict = {}
    designDict = {}
    if paramColn:
      while True:
        line = dIn.readline()
        if line:
          lspl = [x.strip() for x in line.strip().split(',')]
          lspl[0] = [lspl[0]+Postfix[1], lspl[0]+'_R'+lspl[1]]
          lspl[0] = [trackhub.helpers.sanitize(lspl[0][0].replace(".", "_"), strict=False),trackhub.helpers.sanitize(lspl[0][1].replace(".", "_"), strict=False)]
          sampleDesignDict[lspl[0][0]] = {}
          sampleDesignDict[lspl[0][1]] = {}
          for k in paramColn.keys():
            sampleDesignDict[lspl[0][0]][k]=lspl[paramColn[k]]
            sampleDesignDict[lspl[0][1]][k]=lspl[paramColn[k]]
            if k in designDict:
              designDict[k][lspl[paramColn[k]]] = lspl[paramColn[k]]
            else:
              designDict[k] = {lspl[paramColn[k]]:lspl[paramColn[k]]}
        else:
          break

    
    dIn.close()
    
    fileList = []
    fin = open(ListFile,'r')
    while True:
        line = fin.readline()
        if line:
            ifile = line.strip()
            colour = ""
            if sampleDesignDict:
                kfile = trackhub.helpers.sanitize(os.path.splitext(os.path.basename(ifile))[0].replace(".", "_"), strict=False)
                if kfile in sampleDesignDict:
                  if "color" in sampleDesignDict[kfile]:
                    h = sampleDesignDict[kfile]["color"].lstrip('#')
                    colour = ','.join(str(x) for x in tuple(int(h[i:i+2], 16) for i in (0, 2, 4)))
            if len(colour.strip()) == 0:
              colour = '0,0,178'
            fileList.append((PathPrefix.strip()+ifile,colour))
        else:
            break
            fin.close()

    # Initialize the components of a track hub, already connected together
    hub, genomes_file, genome, trackdb = trackhub.default_hub(
        hub_name="RNISRS_hub",
        short_label='Regeneromics Shared Resource hub',
        long_label='Regeneration Next Initiative Regeneromics Shared Resource hub',
        genome=Genome,
        email=EMAIL)
    
    # create compositeTracks
    if sampleDesignDict:
      composite = trackhub.CompositeTrack(
        name = 'composite',
        short_label='singlal'
      )
      # Add those subgroups to the composite track
      subgroups = []
      for k in designDict.keys():
        if k!='color':
          subg = trackhub.SubGroupDefinition(
            name=k,
            label=k,
            mapping=designDict[k]
          )
          subgroups.append(subg)
      
      composite.add_subgroups(subgroups)
      
      # Add the composite track to the trackDb
      trackdb.add_tracks(composite)
      signal_view = trackhub.ViewTrack(
          name='signalviewtrack',
          view='signal',
          short_label='Signal')
      composite.add_view(signal_view)
      regions_view = trackhub.ViewTrack(
          name='regionsviewtrack',
          view='regions',
          short_label='Regions')
      composite.add_view(regions_view)
    
    for ifile,color in fileList:
        extension = os.path.splitext(ifile)[1].replace(".", "").lower()
        filename = trackhub.helpers.sanitize(os.path.splitext(os.path.basename(ifile))[0].replace(".", "_"), strict=False)
        if extension in ['bed','broadpeak','narrowpeak']:
          pass
        elif extension in TrackType.keys():
          if sampleDesignDict:
            track = trackhub.Track(
              name=filename,
              source=ifile,
              color=color,
              visibility=Visibility[extension],
              tracktype=TrackType[extension],
              subgroups=sampleDesignDict[filename],
              autoScale='on')
            signal_view.add_tracks(track)
          else:
            track = trackhub.Track(
              name=filename,
              source=ifile,
              color=color,
              visibility=Visibility[extension],
              tracktype=TrackType[extension],
              autoScale='on')
            trackdb.add_tracks(track)
          linkname=os.path.join(OutFolder, Genome, filename+"."+TrackType[extension])
          makedir(os.path.join(OutFolder, Genome))
          os.symlink(ifile, linkname)
        else:
          pass
    
    hub.render(staging=OutFolder)

############################################
############################################
## RUN FUNCTION
############################################
############################################
def main(args=None):
    args = parse_args(args)
    create_trackhub(
      OutFolder=args.OUTPUT_FOLDER,
      ListFile=args.LIST_FILE,
      Genome=args.GENOME,
      EMAIL=args.EMAIL,
      DesignFile=args.DESIGN_FILE,
      Postfix=args.POSTFIX.split('__'),
      PathPrefix=args.PATH_PREFIX)

if __name__ == "__main__":
    sys.exit(main())


