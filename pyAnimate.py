import os
from subprocess import run
import argparse


parser = argparse.ArgumentParser(description='Create an MP4 from individual PNG files using FFMPEG')
parser.add_argument("-i", "--inDir", help="path to directory containing source PNG files",default='')
parser.add_argument("-f", "--frmRate", help="output frame rate in fps",default='10')
parser.add_argument("-g", "--inGlob", help="glob pattern to identify input files",default='*.png')
parser.add_argument("-c", "--crf", help="constant rate factor, 0-51. 0=lossless, 51=worst quality",default='25')
parser.add_argument("-o", "--outFn", help="output filename",default='out.mp4')

args = parser.parse_args()

inDir = args.inDir
if inDir == '':
    inDir = os.getcwd()
os.chdir(inDir)

frmRate  = args.frmRate
inGlob  = args.inGlob
crf  = args.crf
outFn  = args.outFn

run(['ffmpeg -r ' + frmRate + ' -f image2 -s 1920x1080 -pattern_type glob -i "' + inGlob + '" -vcodec libx264 -crf ' + crf + ' -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" ' + outFn], shell=True, check=True) 