import re
import csv
import os
import argparse
from functionsclinvar import readresultsfkrp, percent_pathogenic, percent_vus

parser = argparse.ArgumentParser(description="sets file name")
parser.add_argument("startfile", help="Enter name of data file")
parser.add_argument("edits", help="Enter name of edits file")
parser.add_argument("destination", help="Enter name of destination file")
args = parser.parse_args()

if os.path.exists(args.startfile):

    with open(args.startfile) as tsv:
        reader = csv.DictReader(tsv, delimiter="\t")

        errors, results = readresultsfkrp(reader, 'FKRP', '077277')
        numerator_path, denom_path, percentage_path = percent_pathogenic(results, 334, 374)
        numerator_vus, denom_vus, percentage_vus = percent_vus(results, 334, 374)

        print 'For pathogenic:'
        print 'On TH: ' + str(numerator_path)
        print 'Total: ' + str(denom_path)
        print str(numerator_path) + '/' + str(denom_path) + ' = ' + str(percentage_path) + '%'

        print 'For VUS:'
        print 'On TH: ' + str(numerator_vus)
        print 'Total: ' + str(denom_vus)
        print str(numerator_vus) + '/' + str(denom_vus) + ' = ' + str(percentage_vus) + '%'

        with open(args.destination, 'w') as f:
            for error in errors:
                f.write('<!--' + error + '-->\n')

            f.write('<svg height="100" width="1500" xmlns="http://www.w3.org/2000/svg">\n\n')

            height = 50
            length = 495 * 3
            # LicD
            f.write('\n<!-- LicD -->\n')
            f.write(
                '\t<polygon points="1003,' + str(height - 6) + ' 1123,' + str(height - 6) + ' 1123,'
                + str(height + 6) + ' 1003,' + str(height + 6) + '" '
                'style = "fill: #9bfffb; fill-opacity: 0.6; stroke: #22e2db; stroke-width:1" />\n\n'
            )
            # line and labels
            f.write('\n<!-- #LINE + LABELS -->\n')
            f.write(
                '\t<line x1="0" y1="' + str(height) + '" x2="'+ str(length) + '" y2="' + str(height) + '" '
                'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
            )
            f.write(
                '\t<line x1="1" y1="' + str(height - 10) + '" x2="1" y2="' + str(height + 10) + '" '
                'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
            )
            f.write(
                '\t<line x1="' + str(length+2) + '" y1="' + str(height - 10) + '" x2="'+ str(length+2) + '" y2="' + str(base + 10) + '" '
                'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
            )
            # ticks (every 100 AAs)
            f.write('\n<!-- TICKS ON #LINE -->\n')
            position = 1
            while position <= length:
                f.write(
                    '\t<line x1="' + str(position) + '" y1="' + str(height - 5) +
                    '" x2="' + str(position) + '" y2="' + str(height + 5) + '" '
                    'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
                )
                position += 300


            # triangle markers

            f.write('\n<!-- MARKERS -->\n')
            for sublist in results:
                num = sublist[0]
                start = str(num * 3 + 1 - 7)
                mid = str(num * 3 + 1)
                end = str(num * 3 + 1 + 7)
                fill = ""
                stroke = ""

                if sublist[1] == 0:
                    print sublist
                    fill = '#bababa'  # gray
                    stroke = '#545454'

                    f.write(
                        '\t<polygon points="' + start + ',20 ' + end + ',20 ' + mid + ',35" '
                                                                                        'style = "fill: ' + fill + '; stroke: ' + stroke + '; stroke-width:1" />\n'
                    )
                    f.write(
                        '\t<line x1="' + mid + '" y1="34" x2="' + mid + '" y2="49" '
                                                                         'style="stroke: ' + stroke + '; stroke-width:2" />\n'
                    )

            f.write('</svg>\n')
else:
    print "error"

# checks:
# command line stuff
# protein is the Right One tm
# spot check
# stick one at the end
# stick one outside bounds
# add up
# make stuff functions
