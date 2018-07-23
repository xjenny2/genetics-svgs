import re
import csv
import os
import argparse
from functionsclinvar import readresultsascending, percent_pathogenic, percent_vus

parser = argparse.ArgumentParser(description="sets file name")
parser.add_argument("startfile", help="Enter name of data file")
parser.add_argument("destination", help="Enter name of destination file")
args = parser.parse_args()

if os.path.exists(args.startfile):

    with open(args.startfile) as tsv:
        reader = csv.DictReader(tsv, delimiter="\t")

        errors, results = readresultsascending(reader, 'COL6A1', '001839')
        numerator_path, denom_path, percentage_path = percent_pathogenic(results, 253, 564)
        numerator_vus, denom_vus, percentage_vus = percent_vus(results, 253, 564)

        print errors
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

            f.write('<svg height="100" width="1600" xmlns="http://www.w3.org/2000/svg">\n\n')
            # triple helical zoom
            f.write('\n<!-- #LINE + LABELS FOR ZOOM-->\n')
            th_start_raw = 253 * 5 + 2
            th_start = 0
            th_end = (564 * 5 + 2) - th_start_raw
            f.write(
                '\t<polygon points="' + str(th_start) + ',38 ' + str(th_end) + ',38 ' + str(th_end) + ',62 ' + str(
                    th_start) + ',62" '
                             'style = "fill: #9bfffb; fill-opacity: 0.6; stroke: #22e2db; stroke-width:2" />\n\n'
            )
            f.write(
                '\t<line x1="' + str(th_start) + '" y1="50" x2="' + str(th_end) + '" y2="50" '
                                                                            'style="stroke: #000000; stroke-width:3; stroke-linecap: round" />\n'
            )
            position = 235
            while th_start <= position <= th_end:
                f.write(
                    '\t<line x1="' + str(position) + '" y1="40" x2="' + str(position) + '" y2="60" '
                                                                        'style="stroke: #000000; stroke-width:3; stroke-linecap: round" />\n'
                )
                position += 500
            for sublist in results:
                num = sublist[0]
                fill = ""
                stroke = ""
                if 253 <= num <= 564 and sublist[1] == 0:
                    print sublist
                    start = str((num * 5 + 2) - 7 - th_start_raw)
                    mid = str((num * 5 + 2) - th_start_raw)
                    end = str((num * 5 + 2) + 7 - th_start_raw)
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
