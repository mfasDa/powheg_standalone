#! /usr/bin/env python3
import argparse

def create_powheginput(origfile, mufact, muren, tag):
    with open(origfile, 'r') as reader and open("powheg.input" as writer):
        for line in reader:
            if "renscfact" in line:
                writer.write("renscfact %d\n" %muren)
                continue
            if "facscfact" in line:
                writer.write("facscfact %d\n" %mufact)
                continue
            if "lhrwgt_id" in line:
                writer.write("lhrwgt_id %s" %tag)
                continue
            writer.write(line)
        reader.close()
        writer.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser("createScaleVariation.py", help="Tool creating POWHEG input for scale variations")
    parser.add_argument("-i", "--inputfile", type=str, required=True, help="Reference configuration")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Reweight tag")
    parser.add_argument("-r", "--renscale", type=float, default=1.0, help="Renormalisation scale")
    parser.add_argument("-f", "--factscale", type=float, default=1.0, help="Factorization scale")
    args = parser.parse_args()
    create_powheginput(args.inputfile, args.factscale, args.renscale, args.tag)