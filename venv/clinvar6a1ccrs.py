import re
import csv
import os
import argparse
import functionsclinvar as fc
import shapes
parser = argparse.ArgumentParser(description="sets file name")
parser.add_argument("startfile", help="Enter name of data file")
parser.add_argument("edits", help="Enter name of edits file")
parser.add_argument("destination", help="Enter name of destination file")
args = parser.parse_args()

if os.path.exists(args.startfile) and os.path.exists(args.edits) and os.path.exists(args.destination):
    print "All files exist--running program"


    with open(args.startfile) as tsv:
        reader = csv.DictReader(tsv, delimiter="\t")

        resultsraw, errors = fc.readresults(reader, 'COL6A1', '001839')
        results, errors = fc.checkascending(resultsraw, errors)
        numerator_path, denom_path, percentage_path = fc.percent_pathogenic(results, 253, 564)
        numerator_vus, denom_vus, percentage_vus = fc.percent_vus(results, 253, 564)

        print errors
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
        edits_list = fc.edits(editsfile, '001839', 6)
    with open("/Users/jennyxu/Desktop/genetics/col6a1_ccrs.tsv") as tsv:
        ccrsfile = csv.reader(tsv, delimiter="\t")
        ccrslist = fc.ccrs(ccrsfile)

        with open(args.destination, 'w') as f:
            for error in errors:
                f.write('<!--' + error + '-->\n')

            f.write('<svg height="500" width="3100" xmlns="http://www.w3.org/2000/svg">\n\n')
            height = 100
            length = 1028 * 3
            # triple helical region
            f.write('\n<!-- TH -->\n')
            f.write(shapes.bluerect(height, 758, 1691))

            # VWA
            f.write('\n<!-- VWA -->\n')
            f.write(shapes.greenrect(height, 110, 683))

            # VWA
            f.write('\n<!-- VWA -->\n')
            f.write(shapes.greenrect(height, 1844, 2363))

            # VWA
            f.write('\n<!-- VWA -->\n')
            f.write(shapes.greenrect(height, 2586, 3056))

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
                '\t<line x1="' + str(length+2) + '" y1="' + str(height - 10) + '" x2="'+ str(length+2) + '" y2="'
                + str(height + 10) + '" style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
            )

            # ticks (every 100 AAs)
            f.write('\n<!-- TICKS ON #LINE -->\n')
            position = 299
            while position <= length:
                f.write(
                    '\t<line x1="' + str(position) + '" y1="' + str(height - 5) +
                    '" x2="' + str(position) + '" y2="' + str(height + 5) + '" '
                    'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
                )
                position += 300

            # triangle markers
            f.write('\n<!-- MARKERS -->\n')

            list_location = []
            list_location_path = []

            for sublist in results:
                overlaps = 0
                overlaps_path = 0

                num = sublist[0]
                overlaps = fc.checkoverlaps(num, list_location, overlaps)
                list_location.append(num)
                if sublist[1] != 0:
                    overlaps_path = fc.checkoverlaps(num, list_location_path, overlaps_path)
                    list_location_path.append(num)
                fill = ""
                stroke = ""
                if sublist[1] == 2 or sublist[1] == 1:
                    fill = '#FF6D8F'  # red
                    stroke = '#960c2c'

                    if sublist[2] == 2:  # nonsense
                        fill = '#FFFF8E'  # yellow
                        stroke = '#ffaa00'

                    elif sublist[2] == 1:  # gly
                        fill = '#91bdff'  # blue
                        stroke = '#004ec4'

                    if overlaps_path == 0:
                        f.write(shapes.triangle(height, num, fill, stroke))
                        f.write(shapes.line(height, num, stroke))
                    else:
                        line2 = height - 30 - (overlaps_path * 9)
                        f.write(shapes.overlap_triangle(line2, num, fill, stroke))
                        
            for sublist in edits_list:  # editable bases
                mid1 = str(int(sublist) * 3 - 1)
                f.write(shapes.circle(height - 35, mid1))

            # rectangle for ccr %s
            f.write('<polygon points="0,170 ' + str(length+2) + ',170 ' + str(length+2) + ',200 0,200" />\n')

            for sublist in ccrslist:  # [position, #]
                ccrspos = str(sublist[0] * 3 - 1) # scale
                percent = float(sublist[1]) / float(100)
                if sublist[1] == -1:  # plots white if no data
                    f.write(
                        '\t<line x1="' + ccrspos + '" y1="172" x2="' + ccrspos + '" y2="198" '
                        'style="stroke: white; stroke-width: 4" />\n'
                    )
                else: # plots w/ color scale
                    f.write(
                    '\t<line x1="' + ccrspos + '" y1="172" x2="' + ccrspos + '" y2="198" '
                    'style="stroke: rgb(' + str(int((percent * 255))) + ',0,'
                    + str(int((255 - (percent * 255)))) + '); stroke-width: 4" />\n'
                    )
            f.write('</svg>\n')
else:
    print "error"
