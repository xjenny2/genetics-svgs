import re
import csv
import os
import argparse
from functionsclinvar import readresultsispd, percent_pathogenic, percent_vus, editsispd, ccrs

parser = argparse.ArgumentParser(description="sets file name")
parser.add_argument("startfile", help="Enter name of data file")
parser.add_argument("edits", help="Enter name of edits file")
parser.add_argument("destination", help="Enter name of destination file")
args = parser.parse_args()

if os.path.exists(args.startfile):

    with open(args.startfile) as tsv:
        reader = csv.DictReader(tsv, delimiter="\t")

        errors, results = readresultsispd(reader, 'ISPD', '001094896')
        numerator_path, denom_path, percentage_path = percent_pathogenic(results, 47, 279)
        numerator_vus, denom_vus, percentage_vus = percent_vus(results, 47, 279)

        print 'For pathogenic:'
        print 'On TH: ' + str(numerator_path)
        print 'Total: ' + str(denom_path)
        print str(numerator_path) + '/' + str(denom_path) + ' = ' + str(percentage_path) + '%'

        print 'For VUS:'
        print 'On TH: ' + str(numerator_vus)
        print 'Total: ' + str(denom_vus)
        print str(numerator_vus) + '/' + str(denom_vus) + ' = ' + str(percentage_vus) + '%'

    with open(args.edits) as tsv:
        editsfile = csv.reader(tsv, delimiter="\t")
        edits_list = editsispd(editsfile, '001094896')

    with open("/Users/jennyxu/Desktop/genetics/ISPD_ccrs") as tsv:
        ccrsfile = csv.reader(tsv, delimiter="\t")
        ccrslist = ccrs(ccrsfile)

        with open(args.destination, 'w') as f:

            for error in errors:
                f.write('<!--' + error + '-->\n')

            f.write('<svg height="200" width="1500" xmlns="http://www.w3.org/2000/svg">\n\n')

            base = 50
            length = 451 * 3
            # IspD
            f.write('\n<!-- IspD -->\n')
            f.write(
                '\t<polygon points="142,' + str(base - 6) + ' 838,' + str(base - 6) + ' 838,'
                + str(base + 6) + ' 142,' + str(base + 6) + '" '
                'style = "fill: #9bfffb; fill-opacity: 0.6; stroke: #22e2db; stroke-width:1" />\n\n'
            )
            # line and labels
            f.write('\n<!-- #LINE + LABELS -->\n')
            f.write(
                '\t<line x1="0" y1="' + str(base) + '" x2="'+ str(length) + '" y2="' + str(base) + '" '
                'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
            )
            f.write(
                '\t<line x1="1" y1="' + str(base - 10) + '" x2="1" y2="' + str(base + 10) + '" '
                'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
            )
            f.write(
                '\t<line x1="'+ str(length + 2) + '" y1="' + str(base - 10) + '" x2="'+ str(length + 2) + '" y2="' + str(base + 10) + '" '
                'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
            )
            # ticks (every 100 AAs)
            f.write('\n<!-- TICKS ON #LINE -->\n')
            position = 1
            while position <= length:
                f.write(
                    '\t<line x1="' + str(position) + '" y1="' + str(base - 5) +
                    '" x2="' + str(position) + '" y2="' + str(base + 5) + '" '
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

                if sublist[1] == 2 or sublist[1] == 1:
                    print sublist
                    fill = '#ffa3aa'  # pink
                    stroke = '#ff5e6b'
                    f.write(
                        '\t<polygon points="' + start + ',20 ' + end + ',20 ' + mid + ',35" '
                        'style = "fill: ' + fill + '; stroke: ' + stroke + '; stroke-width:1" />\n'
                    )
                    f.write(
                        '\t<line x1="' + mid + '" y1="34" x2="' + mid + '" y2="49" '
                        'style="stroke: ' + stroke + '; stroke-width:2" />\n'
                    )

                if (sublist[1] == 1 or sublist[1] == 2) and sublist[2] == 1:
                    fill = '#ffde72'  # yellow
                    stroke = '#ffa716'
                    f.write(
                        '\t<polygon points="' + start + ',20 ' + end + ',20 ' + mid + ',35" '
                                                                                        'style = "fill: ' + fill + '; stroke: ' + stroke + '; stroke-width:1" />\n'
                    )
                    f.write(
                        '\t<line x1="' + mid + '" y1="34" x2="' + mid + '" y2="49" '
                                                                         'style="stroke: ' + stroke + '; stroke-width:2" />\n'
                    )
            for sublist in edits_list:
                editpos = str(int(sublist) * 3 + 1)
                f.write(
                    '<circle cx="' + pos + '" cy="12" r="4"  stroke-width=".5" />'
                )

            f.write('<polygon points="0,70 ' + str(length + 2) + ',70 ' + str(length + 2) + ',100 0,100" />\n')
            for sublist in ccrslist:
                print sublist
                ccrspos = str(sublist[0] * 3 + 1)
                percent = float(sublist[1]) / float(100)
                if sublist[1] == -1:
                    print sublist
                    f.write(
                        '\t<line x1="' + ccrspos + '" y1="72" x2="' + ccrspos + '" y2="98" '
                                                                                'style="stroke: black; stroke-width: 3" />\n'
                    )
                else:
                    f.write(
                    '\t<line x1="' + ccrspos + '" y1="72" x2="' + ccrspos + '" y2="98" '
                    'style="stroke: rgb(' + str(int((percent * 255))) + ',0,' + str(int((255 - (percent * 255)))) + '); stroke-width: 3" />\n'
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
