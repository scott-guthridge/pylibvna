#!/usr/bin/python3
import argparse
from libvna.data import NPData, FileType

usage = """
    %s [-f format] input-file output-file
    where format is a comma-separated list of:
      s[ri|ma|dB]  scattering parameters
      t[ri|ma|dB]  scattering-transfer parameters
      u[ri|ma|dB]  inverse scattering-transfer parameters
      z[ri|ma]     impedance parameters
      y[ri|ma]     admittance parameters
      h[ri|ma]     hybrid parameters
      g[ri|ma]     inverse-hybrid parameters
      a[ri|ma]     ABCD parameters
      b[ri|ma]     inverse ABCD parameters
      Zin[ri|ma]   input impedances
      PRC          Zin as parallel RC
      PRL          Zin as parallel RL
      SRC          Zin as series RC
      SRL          Zin as series RL
      IL           insertion loss
      RL           return loss
      VSWR         voltage standing wave ratio

    Coordinates
      ri  real, imaginary
      ma  magnitude, angle
      dB  decibels, angle

    Specifiers are case-insensitive."

    File Types
      Touchstone v1:          .s1p, .s2p, .s3p, .s4p
      Touchstone v2:          .ts
      Network Parameter Data: .npd
"""

parser = argparse.ArgumentParser(usage=usage)
parser.add_argument("-f", "--format")
parser.add_argument("input_file", metavar="input-file")
parser.add_argument("output_file", metavar="output-file")
args = parser.parse_args()

npd = NPData()
npd.load(args.input_file)
npd.filetype = FileType.AUTO
if args.format:
    npd.format = args.format
npd.save(args.output_file)
