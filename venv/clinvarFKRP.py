import re
import csv
import os
import argparse
from functionsclinvar import readresultsfkrp, percent_pathogenic, percent_vus

parser = argparse.ArgumentParser(description="sets file name")
parser.add_argument("startfile", help="Enter name of data file")
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
            f.write('<!DOCTYPE html>\n<html>\n')
            f.write('<head></head>\n')
            f.write('<body>\n')

            for error in errors:
                f.write('<!--' + error + '-->\n')

            f.write('<svg height="500" width="2100">\n\n')

            # box for key
            f.write('\n<!-- KEY -->\n')
            f.write(
                '\t<rect x="1" y="1" width="400" height="160" '
                'style="fill: white; stroke-width:1; stroke:rgb(0,0,0)" />\n\n'
            )
            f.write(
                '\t\t<text text-anchor="middle" x="150" y="35" '
                'style = "font-size: 30; font-weight: bold; text-decoration: underline;">FKRP (clinvar)</text>\n'
            )
            f.write(
                '\t\t<text text-anchor="start" x="20" y="65" '
                'style = "font-size: 10">Key:</text>\n\n'
            )
            # key PINK
            f.write(
                '\t\t<rect x="20" y="80" width="10" height="10" '
                'style="fill: FE018C; stroke-width:1; stroke: #60002c" />\n'
            )
            f.write(
                '\t\t<text text-anchor="start" x="40" y="90" '
                'style = "font-size: 10; font-weight: bold; fill: #FE018C">Pathogenic</text>\n\n'
            )
            # key YELLOW
            f.write(
                '\t\t<rect x="20" y="100" width="10" height="10" '
                'style="fill: #fff600; stroke-width:1; stroke: #593e00" />\n'
            )
            f.write(
                '\t\t<text text-anchor="start" x="40" y="110" '
                'style = "font-size: 10; font-weight: bold; fill: #f28d00">Likely Pathogenic</text>\n\n'
            )
            # key GREEN
            f.write(
                '\t\t<rect x="20" y="120" width="10" height="10" '
                'style="fill: #00E408; stroke-width:1; stroke: #005908" />\n'
            )
            f.write(
                '\t\t<text text-anchor="start" x="40" y="130" '
                'style = "font-size: 10; font-weight: bold; fill: #00E408">Uncertain</text>\n\n'
            )
            # key LicD
            f.write(
                '\t\t<rect x="20" y="140" width="10" height="10" '
                'style="fill: #00ABFF; stroke-width:1; stroke: #002f82" />\n'
            )
            f.write(
                '\t\t<text text-anchor="start" x="40" y="150" '
                'style = "font-size: 10; font-weight: bold; fill: #00ABFF">LicD</text>\n\n'
            )
            # key Ter
            f.write(
                '\t\t<rect x="250" y="80" width="10" height="10" '
                'style="fill: #8d02ff; stroke-width:1; stroke: #23004c" />\n'
            )
            f.write(
                '\t\t<text text-anchor="start" x="270" y="90" '
                'style = "font-size: 10; font-weight: bold; fill: #8d02ff">Ter</text>\n\n'
            )


            n = 0
            for n in range(0, 2):
                base = 200 + 100 * n
                # LicD
                f.write('\n<!-- LicD -->\n')
                f.write(
                    '\t<polygon points="668,' + str(base - 6) + ' 748,' + str(base - 6) + ' 748,'
                    + str(base + 6) + ' 668,' + str(base + 6) + '" '
                    'style = "fill: #00ABFF; fill-opacity: 0.6; stroke: #002f82; stroke-width:1" />\n\n'
                )
                # line and labels
                f.write('\n<!-- #LINE + LABELS -->\n')
                f.write(
                    '\t<line x1="0" y1="' + str(base) + '" x2="990" y2="' + str(base) + '" '
                    'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
                )
                f.write(
                    '\t<line x1="1" y1="' + str(base - 10) + '" x2="1" y2="' + str(base + 10) + '" '
                    'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
                )
                f.write(
                    '<text text-anchor="start" x="0" y="' + str(base + 30) + '" '
                    'style = "font-size: 12px;">0</text>'
                )
                f.write(
                    '\t<line x1="990" y1="' + str(base - 10) + '" x2="990" y2="' + str(base + 10) + '" '
                    'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
                )
                f.write(
                    '\t<text text-anchor="end" x="1000" y="' + str(base + 30) + '" '
                    'style = "font-size: 12px;">990</text>\n\n'
                )
                # ticks (every 100 AAs)
                f.write('\n<!-- TICKS ON #LINE -->\n')
                position = 0
                while position <= 990:
                    f.write(
                        '\t<line x1=' + str(position) + ' y1="' + str(base - 5) +
                        '" x2=' + str(position) + ' y2="' + str(base + 5) + '" '
                        'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
                    )
                    position += 200
                n += 1

            # triangle markers
            f.write('\n<!-- MARKERS -->\n')
            for sublist in results:
                num = sublist[0]
                start = str(num * 2 - 7)
                mid = str(num * 2)
                end = str(num * 2 + 7)
                fill = ""
                stroke = ""

                if sublist[1] == 2:
                    fill = '#FE018C'  # pink
                    stroke = '#60002c'
                    f.write(
                        '\t<polygon points="' + start + ',270 ' + end + ',270 ' + mid + ',285" '
                        'style = "fill: ' + fill + '; stroke: ' + stroke + '; stroke-width:1" />\n'
                    )
                    f.write(
                        '\t<line x1="' + mid + '" y1="284" x2="' + mid + '" y2="299" '
                        'style="stroke: ' + stroke + '; stroke-width:2" />\n'
                    )

                elif sublist[1] == 1:
                    fill = '#fff600'  # yellow
                    stroke = '#593e00'
                    f.write(
                        '\t<polygon points="' + start + ',270 ' + end + ',270 ' + mid + ',285" '
                        'style = "fill: ' + fill + '; stroke: ' + stroke + '; stroke-width:1" />\n'
                    )
                    f.write(
                        '\t<line x1="' + mid + '" y1="284" x2="' + mid + '" y2="299" '
                        'style="stroke: ' + stroke + '; stroke-width:2" />\n'
                    )

                elif sublist[1] == 0:
                    fill = '#00E408'  # green
                    stroke = '#005908'

                f.write(
                    '\t<polygon points="' + start + ',170 ' + end + ',170 ' + mid + ',185" '
                    'style = "fill: ' + fill + '; stroke: ' + stroke + '; stroke-width:1" />\n'
                )
                f.write(
                    '\t<line x1="' + mid + '" y1="184" x2="' + mid + '" y2="199" '
                    'style="stroke: ' + stroke + '; stroke-width:2" />\n'
                )
                if (sublist[1] == 1 or sublist[1] == 2) and sublist[2] == 1:
                    fill = '#8d02ff'
                    stroke = '#23004c'
                    f.write(
                        '\t<polygon points="' + start + ',270 ' + end + ',270 ' + mid + ',285" '
                                                                                        'style = "fill: ' + fill + '; stroke: ' + stroke + '; stroke-width:1" />\n'
                    )
                    f.write(
                        '\t<line x1="' + mid + '" y1="284" x2="' + mid + '" y2="299" '
                                                                         'style="stroke: ' + stroke + '; stroke-width:2" />\n'
                    )

            f.write('</svg>\n')
            f.write('</body>\n')
            f.write('</html>')
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
