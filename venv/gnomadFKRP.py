import re
import csv
from functionsgnomad import readresults, checkascending
with open("/Users/jennyxu/Desktop/genetics/gnomad_FKRP") as csvDataFile:
    reader = csv.DictReader(csvDataFile)  # sets up file as dictionary with first line as header
    raw_results = readresults(reader)
    results, errors = checkascending(raw_results)
    print results
    print errors


    with open('/Users/jennyxu/Desktop/genetics/gnomad_FKRP_svg.html', 'w') as f:
        f.write('<!DOCTYPE html>\n<html>\n')
        f.write('<head>\n</head>\n')
        f.write('<body>\n')
        for error in errors:
            f.write('<!--' + error + '-->\n')
        print errors
        f.write('<svg height="1000" width="2100">\n\n')
        # box for key
        f.write('\n<!-- KEY -->\n')
        f.write(
            '\t<rect x="1" y="1" width="560" height="160" '
            'style="fill: white; stroke-width:1; stroke:rgb(0,0,0)" />\n'
        )
        f.write(
            '\t\t<text text-anchor="middle" x="275" y="35" '
            'style = "font-size: 30; font-weight: bold; text-decoration: underline;">FKRP (gnomAD)</text>\n'
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
        # key TER
        f.write(
            '\t\t<rect x="295" y="110" width="10" height="10" '
            'style="fill: #8d02ff; stroke-width:1; stroke: #23004c" />\n'
        )
        f.write(
            '\t\t<text text-anchor="start" x="315" y="120" '
            'style = "font-size: 10; font-weight: bold; fill: #8d02ff">Ter</text>\n\n'
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
        # key LicD
        f.write(
            '\t\t<rect x="295" y="140" width="10" height="10" '
            'style="fill: #00ffe9; stroke-width:1; stroke: #3ca390" />\n'
        )
        f.write(
            '\t\t<text text-anchor="start" x="315" y="150" '
            'style = "font-size: 10; font-weight: bold; fill: #00ffe9">LicD</text>\n\n'
        )
        # set up number line + labels + trip helix
        n = 0
        for n in range(0, 5):
            base = 200 + 100 * n
            # Ter
            f.write('\n<!-- TER -->\n')
            f.write(
                '\t<polygon points="668,' + str(base - 6) + ' 748,' + str(base - 6) + ' 748,'
                + str(base + 6) + ' 668,' + str(base + 6) + '" '
                'style = "fill: #00ffe9; fill-opacity: 0.6; stroke: #3ca390; stroke-width:1" />\n\n'
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
            num = int(sublist[0])
            start = str(num*2 - 7)
            mid = str(num*2)
            end = str(num*2 + 7)
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
            f.write(
                '\t<polygon points="' + start + ',180 ' + end + ',180 ' + mid + ',198" '
                'style = "fill: ' + fill + '; stroke: ' + stroke + '; stroke-width:1" />\n'
            )



        f.write('</svg>\n')
        f.write('</body>\n')
        f.write('</html>')

# checks:
# spot checked code
# NEED TO CHECK MATH FOR ZOOMED IN FIGURE
