import re
import csv

with open("/Users/jennyxu/Desktop/genetics/gnomad_col6a3") as csvDataFile:
    reader = csv.DictReader(csvDataFile)  # sets up file as dictionary with first line as header
    results = []
    just_location = []
    errors = []
    pattern = re.compile(
        r'(?<=p\.[A-Z][a-z]{2})\d+(?=[A-Z][a-z]{2}$)'  # searches for p.XXX[Location]XXX"
    )
    intpattern = re.compile(r'[\d.]+(?=e-)')  # finds base in exponent
    intexp = re.compile(r'(?<=e)[-\d]+')  # finds exponent
    glypattern = re.compile(r'(?<=p.Gly)\d+(?=[A-Z][a-z]+)')  # finds exponent
    for row in reader:
        filter_ex = row["Filters - exomes"]
        filter_gen = row["Filters - genomes"]
        alleleCt = row["Allele Count"]
        alleleFreq = row["Allele Frequency"]
        consequence = row["Protein Consequence"]
        if (filter_ex == "PASS" or filter_gen == "PASS") and pattern.findall(consequence):
            location = str(pattern.findall(consequence))[2:-2]  # remove [''] from string
            location = int(location)
            gly = 0
            if glypattern.findall(consequence):
                gly = 1
            if alleleCt == '1':
                alleleCt = 2
            elif intpattern.findall(alleleFreq):
                freqint = float(str(intpattern.findall(alleleFreq))[2:-2])
                freqexp = float(str(intexp.findall(alleleFreq))[2:-2])
                # alleleFreq = freqint ** (-freqexp)
                alleleCt = freqint * 10 ** freqexp
            else:
                alleleCt = alleleFreq
            alleleCt = float(alleleCt)
            result = [location, alleleCt, gly]
            results.append(result)
    print results
    for sublist in results:
            just_location.append(sublist[0])  # appends first # in each pair to just_location
    for index, item in enumerate(just_location):  # checks if in ascending order, if not print error msg
        if item > just_location[index - 1] and index != 0:
            errors.append("Error at row %d: %d not descending" % (index, item))

    with open('/Users/jennyxu/Desktop/genetics/gnomad_col6A3_svg.html', 'w') as f:
        f.write('<!DOCTYPE html>\n<html>\n')
        f.write('<head>\n</head>\n')
        f.write('<body>\n')
        for error in errors:
            f.write('<!--' + error + '-->\n')
        print errors
        f.write('<svg height="1000" width="3200">\n\n')
        # box for key
        f.write('\n<!-- KEY -->\n')
        f.write(
            '\t<rect x="1" y="1" width="675" height="160" '
            'style="fill: white; stroke-width:1; stroke:rgb(0,0,0)" />\n'
        )
        f.write(
            '\t\t<text text-anchor="middle" x="275" y="35" '
            'style = "font-size: 30; font-weight: bold; text-decoration: underline;">COL6A3 (gnomAD)</text>\n'
        )
        f.write(
            '\t\t<text text-anchor="start" x="20" y="60" '
            'style = "font-size: 10">Key:</text>\n\n'
        )
        # key PINK
        f.write(
            '\t\t<rect x="20" y="80" width="10" height="10" '
            'style="fill: FE018C; stroke-width:1; stroke: #60002c" />\n'
        )
        f.write(
            '\t\t<text text-anchor="start" x="40" y="90" '
            'style = "font-size: 10; font-weight: bold; fill: #FE018C">1 Allele Found</text>\n\n'
        )
        # key YELLOW
        f.write(
            '\t\t<rect x="20" y="110" width="10" height="10" '
            'style="fill: #fff600; stroke-width:1; stroke: #593e00" />\n'
        )
        f.write(
            '\t\t<text text-anchor="start" x="40" y="120" '
            'style = "font-size: 10; font-weight: bold; fill: #f28d00">Freq &le; .25%</text>\n\n'
        )
        # key GREEN
        f.write(
            '\t\t<rect x="20" y="140" width="10" height="10" '
            'style="fill: #00E408; stroke-width:1; stroke: #005908" />\n'
        )
        f.write(
            '\t\t<text text-anchor="start" x="40" y="150" '
            'style = "font-size: 10; font-weight: bold; fill: #00E408">.25% < Freq &le; 1%</text>\n\n'
        )
        # key PURPLE
        f.write(
            '\t\t<rect x="295" y="110" width="10" height="10" '
            'style="fill: #8d02ff; stroke-width:1; stroke: #23004c" />\n'
        )
        f.write(
            '\t\t<text text-anchor="start" x="315" y="120" '
            'style = "font-size: 10; font-weight: bold; fill: #8d02ff">Gly</text>\n\n'
        )
        # key BLUE
        f.write(
            '\t\t<rect x="295" y="80" width="10" height="10" '
            'style="fill: #00ABFF; stroke-width:1; stroke: #002f82" />\n'
        )
        f.write(
            '\t\t<text text-anchor="start" x="315" y="90" '
            'style = "font-size: 10; font-weight: bold; fill: #00ABFF">Freq > 1%</text>\n\n'
        )
        # key TRIPLE HELICAL
        f.write(
            '\t\t<rect x="295" y="140" width="10" height="10" '
            'style="fill: #00ffe9; stroke-width:1; stroke: #3ca390" />\n'
        )
        f.write(
            '\t\t<text text-anchor="start" x="315" y="150" '
            'style = "font-size: 10; font-weight: bold; fill: #00ffe9">Triple Helical Struct.</text>\n\n'
        )
        # key VWA
        f.write(
            '\t\t<rect x="570" y="80" width="10" height="10" '
            'style="fill: #ff9a32; stroke-width:1; stroke: #e57400" />\n'
        )
        f.write(
            '\t\t<text text-anchor="start" x="590" y="90" '
            'style = "font-size: 10; font-weight: bold; fill: #ff9a32">VWA</text>\n\n'
        )
        # set up number line + labels + trip helix
        n = 0
        for n in range(0, 5):
            base = 200 + 100 * n
            f.write('\n<!-- TRIPLE HELICAL REGION -->\n')
            f.write(
                '\t<polygon points="2036,' + str(base - 6) + ' 2096,' + str(base - 6) + ' 2096,'
                + str(base + 6) + ' 2036,' + str(base + 6) + '" '
                'style = "fill: #00ffe9; fill-opacity: 0.6; stroke: #3ca390; stroke-width:1" />\n\n'
            )
            # VWA
            f.write('\n<!-- VWA -->\n')
            f.write(
                '\t<polygon points="39,' + str(base - 6) + ' 212,' + str(base - 6) + ' 212,'
                + str(base + 6) + ' 39,' + str(base + 6) + '" '
                'style = "fill: #ff9a32; fill-opacity: 0.6; stroke: #e57400; stroke-width:1" />\n\n'
            )
            # VWA
            f.write('\n<!-- VWA -->\n')
            f.write(
                '\t<polygon points="242,' + str(base - 6) + ' 414,' + str(base - 6) + ' 414,'
                + str(base + 6) + ' 242,' + str(base + 6) + '" '
                'style = "fill: #ff9a32; fill-opacity: 0.6; stroke: #e57400; stroke-width:1" />\n\n'
            )
            # VWA
            f.write('\n<!-- VWA -->\n')
            f.write(
                '\t<polygon points="445,' + str(base - 6) + ' 616,' + str(base - 6) + ' 616,'
                + str(base + 6) + ' 445,' + str(base + 6) + '" '
                'style = "fill: #ff9a32; fill-opacity: 0.6; stroke: #e57400; stroke-width:1" />\n\n'
            )
            # VWA
            f.write('\n<!-- VWA -->\n')
            f.write(
                '\t<polygon points="639,' + str(base - 6) + ' 811,' + str(base - 6) + ' 811,'
                + str(base + 6) + ' 639,' + str(base + 6) + '" '
                'style = "fill: #ff9a32; fill-opacity: 0.6; stroke: #e57400; stroke-width:1" />\n\n'
            )
            # VWA
            f.write('\n<!-- VWA -->\n')
            f.write(
                '\t<polygon points="837,' + str(base - 6) + ' 1008,' + str(base - 6) + ' 1008,'
                + str(base + 6) + ' 837,' + str(base + 6) + '" '
                'style = "fill: #ff9a32; fill-opacity: 0.6; stroke: #e57400; stroke-width:1" />\n\n'
            )
            # VWA
            f.write('\n<!-- VWA -->\n')
            f.write(
                '\t<polygon points="1029,' + str(base - 6) + ' 1200,' + str(base - 6) + ' 1200,'
                + str(base + 6) + ' 1029,' + str(base + 6) + '" '
                'style = "fill: #ff9a32; fill-opacity: 0.6; stroke: #e57400; stroke-width:1" />\n\n'
            )
            # VWA
            f.write('\n<!-- VWA -->\n')
            f.write(
                '\t<polygon points="1233,' + str(base - 6) + ' 1403,' + str(base - 6) + ' 1403,'
                + str(base + 6) + ' 1233,' + str(base + 6) + '" '
                'style = "fill: #ff9a32; fill-opacity: 0.6; stroke: #e57400; stroke-width:1" />\n\n'
            )

            # VWA
            f.write('\n<!-- VWA -->\n')
            f.write(
                '\t<polygon points="1436,' + str(base - 6) + ' 1608,' + str(base - 6) + ' 1608,'
                + str(base + 6) + ' 1436,' + str(base + 6) + '" '
                'style = "fill: #ff9a32; fill-opacity: 0.6; stroke: #e57400; stroke-width:1" />\n\n'
            )
            # VWA
            f.write('\n<!-- VWA -->\n')
            f.write(
                '\t<polygon points="1639,' + str(base - 6) + ' 1811,' + str(base - 6) + ' 1811,'
                + str(base + 6) + ' 1639,' + str(base + 6) + '" '
                'style = "fill: #ff9a32; fill-opacity: 0.6; stroke: #e57400; stroke-width:1" />\n\n'
            )
            # VWA
            f.write('\n<!-- VWA -->\n')
            f.write(
                '\t<polygon points="1838,' + str(base - 6) + ' 2003,' + str(base - 6) + ' 2003,'
                + str(base + 6) + ' 1838,' + str(base + 6) + '" '
                'style = "fill: #ff9a32; fill-opacity: 0.6; stroke: #e57400; stroke-width:1" />\n\n'
            )
            # line and labels
            f.write('\n<!-- #LINE + LABELS -->\n')
            f.write(
                '\t<line x1="0" y1="' + str(base) + '" x2="3177" y2="' + str(base) + '" '
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
                '\t<line x1="3177" y1="' + str(base - 10) + '" x2="3177" y2="' + str(base + 10) + '" '
                'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
            )
            f.write(
                '\t<text text-anchor="end" x="3177" y="' + str(base + 30) + '" '
                'style = "font-size: 12px;">3177</text>\n\n'
            )
            # ticks (every 100 AAs)
            f.write('\n<!-- TICKS ON #LINE -->\n')
            position = 0
            while position <= 3177:
                f.write(
                    '\t<line x1=' + str(position) + ' y1="' + str(base - 5) +
                    '" x2=' + str(position) + ' y2="' + str(base + 5) + '" '
                    'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
                )
                position += 100
            n += 1

        # triple helical zoom
        f.write('\n<!-- #LINE + LABELS FOR ZOOM-->\n')
        f.write(
            '\t<polygon points="10,' + str(700 - 12) + ' 730,' + str(700 - 12) + ' 730,'
            + str(700 + 12) + ' 10,' + str(700 + 12) + '" '
            'style = "fill: #00ffe9; fill-opacity: 0.6; stroke: #3ca390; stroke-width:2" />\n\n'
        )
        f.write(
            '\t<line x1="10" y1="' + str(700) + '" x2="730" y2="' + str(700) + '" '
            'style="stroke: #000000; stroke-width:3; stroke-linecap: round" />\n'
        )
        f.write(
            '<text text-anchor="start" x="1" y="' + str(700 + 30) + '" '
            'style = "font-size: 12px;">2036</text>'
        )
        f.write(
            '\t<text text-anchor="middle" x="730" y="' + str(700 + 30) + '" '
            'style = "font-size: 12px;">2096</text>\n\n'
        )

        # triangle markers
        f.write('\n<!-- MARKERS -->\n')
        for sublist in results:
            num = int(sublist[0])
            start = str(num - 7)
            mid = str(num)
            end = str(num + 7)
            fill = ''
            stroke = ''
            if sublist[1] == 2:
                fill = '#FE018C'  # pink
                stroke = '#60002c'

                f.write(
                    '\t<polygon points="' + start + ',580 ' + end + ',580 ' + mid + ',598" '
                    'style = "fill: ' + fill + '; stroke: '
                    + stroke + '; stroke-opacity: 0.8; stroke-width:1" />\n'
                )
            elif sublist[1] <= 0.0025:
                fill = '#fff600'  # yellow
                stroke = '#593e00'

                f.write(
                    '\t<polygon points="' + start + ',480 ' + end + ',480 ' + mid + ',498" '
                    'style = "fill: ' + fill + '; fill-opacity: 0.6; stroke: ' + stroke
                    + '; stroke-opacity: 0.8; stroke-width:1" />\n'
                )
            elif 0.0025 < sublist[1] <= 0.01:
                fill = '#00E408'  # green
                stroke = '#005908'
                f.write(
                    '\t<polygon points="' + start + ',380 ' + end + ',380 ' + mid + ',398" '
                    'style = "fill: ' + fill + '; fill-opacity: 0.6; stroke: ' + stroke
                    + '; stroke-opacity: 0.8; stroke-width:1" />\n'
                )
            elif 0.01 < sublist[1] < 1:
                fill = '#00ABFF'  # blue
                stroke = '#002f82'
                f.write(
                    '\t<polygon points="' + start + ',280 ' + end + ',280 ' + mid + ',298" '
                    'style = "fill: ' + fill + '; fill-opacity: 0.6; stroke: ' + stroke
                    + '; stroke-opacity: 0.8; stroke-width:1" />\n'
                )
            # use lines as markers instead
            # f.write('\t<line x1="' + mid + '" y1="189" x2="' + mid + '" y2="209" style="stroke: ' + fill + '; '
            #         'stroke-width:4" />\n')
            f.write(
                '\t<polygon points="' + start + ',180 ' + end + ',180 ' + mid + ',198" '
                'style = "fill: ' + fill + '; stroke: ' + stroke + '; stroke-width:1" />\n'
            )
            # zoom region markers

            start1 = str(num * 12 - 24422 - 7)
            mid1 = str(num * 12 - 24422)
            end1 = str(num * 12 - 24422 + 7)

            if 2036 <= int(mid) <= 2096:
                f.write(
                    '\t<polygon points="' + start1 + ',670 ' + end1 + ',670 ' + mid1 + ',685" '
                    'style = "fill: ' + fill + '; stroke: ' + stroke + '; stroke-width:1" />\n'
                )
                f.write(
                    '\t<line x1="' + mid1 + '" y1="684" x2="' + mid1 + '" y2="699" '
                    'style="stroke: ' + stroke + '; stroke-width:2" />\n'
                )
                if sublist[2] == 1:
                    fill = '#8d02ff'  # purple
                    stroke = '#23004c'
                    f.write(
                        '<circle cx="' + mid1 + '" cy="665" r="3" stroke="' + stroke + '" stroke-width=".5" fill="' + fill + '" />'
                    )
        f.write('</svg>\n')
        f.write('</body>\n')
        f.write('</html>')

# checks:
# spot checked code
# NEED TO CHECK MATH FOR ZOOMED IN FIGURE
