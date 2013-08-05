import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('files', metavar='files', type=str, nargs='+',
                    help='files to be converted')

args = parser.parse_args()

if args.files:
    files = args.files
    print "Converting "+str(len(files))+" pdf files to png"
else:
    files = os.listdir('.')

def renamer() :
    for filename in files :
        if filename.endswith('pdf'):
            epsfile = filename.replace(".pdf",".eps")
            pngfile = filename.replace(".pdf",".png")
            os.system("pdftops -eps "+filename)
            os.system("convert "+epsfile+" "+pngfile)
            os.remove(epsfile)
            print "converted " + filename + " to " + pngfile
            
# test the function/module
if __name__ == '__main__':
    renamer()
                
