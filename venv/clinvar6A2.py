import re
import csv
import os
import argparse

parser = argparse.ArgumentParser(description="sets file name")
parser.add_argument("startfile", help="Enter name of data file")
parser.add_argument("destination", help="Enter name of destination file")
args = parser.parse_args()
if os.path.exists(args.startfile):
    with open(args.startfile) as tsv:
        reader = csv.DictReader(tsv, delimiter="\t")
        results = []
        just_location = []
        errors = []
        location_pat = re.compile(
            r'(?<=NP_\d{6}\.\d:p\.[A-Z][a-z][a-z])\d+(?=[A-Z][a-z]{2}$)'  # searches for NP_######.#:p.XXX[Location]XXX"
        )
        protein_pat = re.compile(r'(?<=NP_)\d{6}')  # searches for NP_[protein #]
        glypattern = re.compile(r'(?<=p.Gly)\d+(?=[A-Z][a-z]+)')
        for row in reader:
            gene = row["symbol"]
            pathogenic = row["pathogenic"]
            likely_path = row["likely_pathogenic"]
            uncertain = row["uncertain_significance"]
            consequence = row["hgvs_p"]
            key = 0
            gly = 0
            if gene == 'COL6A2' and int(uncertain) >= 1 and location_pat.findall(consequence):
                location = str(location_pat.findall(consequence))[2:-2]
                location = int(location)
                protein = str(protein_pat.findall(consequence))[2:-2]
                if protein == '001840':
                    key = 0
                else:
                    errors.append("Error: protein " + protein + " not protein 001840")
                if glypattern.findall(consequence):
                    gly = 1
                result = [location, key, gly]
                results.append(result)
            elif gene == 'COL6A2' and int(likely_path) >= 1 and location_pat.findall(consequence):
                location = str(location_pat.findall(consequence))[2:-2]
                location = int(location)
                protein = str(protein_pat.findall(consequence))[2:-2]
                if protein == '001840':
                    key = 1
                else:
                    errors.append("Error: protein " + protein + " not protein 001840")
                if glypattern.findall(consequence):
                    gly = 1
                result = [location, key, gly]
                results.append(result)
            elif gene == 'COL6A2' and int(pathogenic) >= 1 and location_pat.findall(consequence):
                location = str(location_pat.findall(consequence))[2:-2]
                location = int(location)
                protein = str(protein_pat.findall(consequence))[2:-2]
                if protein == '001840':
                    key = 2
                else:
                    errors.append("Error: protein " + protein + " not protein 001840")
                if glypattern.findall(consequence):
                    gly = 1
                result = [location, key, gly]
                results.append(result)
        for sublist in results:
            just_location.append(sublist[0])  # appends first # in each pair to just_location
        for index, item in enumerate(just_location):  # checks if in ascending order, if not print error msg
            if item < just_location[index - 1] and index != 0:
                errors.append("Error at row %d: %d not ascending" % (index, item))
        print errors
        print results

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
                '\t<rect x="1" y="1" width="640" height="160" '
                'style="fill: white; stroke-width:1; stroke:rgb(0,0,0)" />\n\n'
            )
            f.write(
                '\t\t<text text-anchor="middle" x="150" y="35" '
                'style = "font-size: 30; font-weight: bold; text-decoration: underline;">COL6A2 (clinvar)</text>\n'
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
            # key TRIPLE HELICAL
            f.write(
                '\t\t<rect x="20" y="140" width="10" height="10" '
                'style="fill: #00ABFF; stroke-width:1; stroke: #002f82" />\n'
            )
            f.write(
                '\t\t<text text-anchor="start" x="40" y="150" '
                'style = "font-size: 10; font-weight: bold; fill: #00ABFF">Triple Helical Structure</text>\n\n'
            )
            # key GLY
            f.write(
                '\t\t<rect x="340" y="80" width="10" height="10" '
                'style="fill: #8d02ff; stroke-width:1; stroke: #23004c" />\n'
            )
            f.write(
                '\t\t<text text-anchor="start" x="360" y="90" '
                'style = "font-size: 10; font-weight: bold; fill: #8d02ff">Gly</text>\n\n'
            )
            # key VWA
            f.write(
                '\t\t<rect x="340" y="100" width="10" height="10" '
                'style="fill: #ff9a32; stroke-width:1; stroke: #e57400" />\n'
            )
            f.write(
                '\t\t<text text-anchor="start" x="360" y="110" '
                'style = "font-size: 10; font-weight: bold; fill: #ff9a32">VWA</text>\n\n'
            )
            n = 0
            for n in range(0, 2):
                base = 200 + 100 * n
                # triple helical region
                f.write('\n<!-- TRIPLE HELICAL REGION -->\n')
                f.write(
                    '\t<polygon points="508,' + str(base - 6) + ' 1180,' + str(base - 6) + ' 1180,'
                    + str(base + 6) + ' 508,' + str(base + 6) + '" '
                    'style = "fill: #00ABFF; fill-opacity: 0.6; stroke: #002f82; stroke-width:1" />\n\n'
                )
                # VWA
                f.write('\n<!-- VWA -->\n')
                f.write(
                    '\t<polygon points="92,' + str(base - 6) + ' 464,' + str(base - 6) + ' 464,'
                    + str(base + 6) + ' 92,' + str(base + 6) + '" '
                    'style = "fill: #ff9a32; fill-opacity: 0.6; stroke: #e57400; stroke-width:1" />\n\n'
                )
                # VWA
                f.write('\n<!-- VWA -->\n')
                f.write(
                    '\t<polygon points="1230,' + str(base - 6) + ' 1596,' + str(base - 6) + ' 1596,'
                    + str(base + 6) + ' 1230,' + str(base + 6) + '" '
                    'style = "fill: #ff9a32; fill-opacity: 0.6; stroke: #e57400; stroke-width:1" />\n\n'
                )
                # VWA
                f.write('\n<!-- VWA -->\n')
                f.write(
                    '\t<polygon points="1666,' + str(base - 6) + ' 2020,' + str(base - 6) + ' 2020,'
                    + str(base + 6) + ' 1666,' + str(base + 6) + '" '
                    'style = "fill: #ff9a32; fill-opacity: 0.6; stroke: #e57400; stroke-width:1" />\n\n'
                )
                # line and labels
                f.write('\n<!-- #LINE + LABELS -->\n')
                f.write(
                    '\t<line x1="0" y1="' + str(base) + '" x2="2038" y2="' + str(base) + '" '
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
                    '\t<line x1="2038" y1="' + str(base - 10) + '" x2="2038" y2="' + str(base + 10) + '" '
                    'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
                )
                f.write(
                    '\t<text text-anchor="end" x="2062" y="' + str(base + 30) + '" '
                    'style = "font-size: 12px;">1019</text>\n\n'
                )
                # ticks (every 100 AAs)
                f.write('\n<!-- TICKS ON #LINE -->\n')
                position = 0
                while position <= 2038:
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
                    stroke = '#23004c' #purple
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
