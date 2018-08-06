import csv
import re
import pysam
import requests
import shapes
import functionsgnomad as fg

results = []

chrom = raw_input("Chromosome number: ")
startpos = raw_input("Genomic Start Pos: ")
endpos = raw_input("Genomic End Pos: ")

proteinid = raw_input('Enter Protein ID: ')
gene_name = fg.genename(proteinid)


# domains
repeats = raw_input('How many kinds of domains?: ')
while repeats.isdigit() is False:
    print "Error: not a number"
    repeats = raw_input('How many kinds of domains? ')
# prints out info for check
domainresults, protlength_raw = fg.domains(repeats, proteinid)
print 'Gene: ' + gene_name
print 'Length: ' + str(protlength_raw) + ' amino acids'

# 47401651
# 47424964
final_results = fg.readresultsgeneral(chrom, startpos, endpos)
results, errors = fg.checkascending(final_results)
if len(results) == 0:
    print "Error: no results found"
    raise SystemExit

with open('/Users/jennyxu/Desktop/genetics/gnomad_col6A1_rare.svg', 'w') as f:
    f.write('<svg height="1000" width="' + str(protlength_raw * 2 + 30) + '" xmlns="http://www.w3.org/2000/svg">\n\n')
    # box for key
    f.write('\n<!-- KEY -->\n')
    f.write(
        '\t<rect x="1" y="1" width="560" height="160" '
        'style="fill: white; stroke-width:1; stroke:rgb(0,0,0)" />\n'
    )
    f.write(
        '\t\t<text text-anchor="middle" x="275" y="35" '
        'style = "font-size: 30; font-weight: bold; text-decoration: underline;">' + gene_name + ' (gnomAD)</text>\n'
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
        'style = "font-size: 10; font-weight: bold; fill: #f28d00">Freq &#60; .25%</text>\n\n'
    )
    # key GREEN
    f.write(
        '\t\t<rect x="20" y="140" width="10" height="10" '
        'style="fill: #00E408; stroke-width:1; stroke: #005908" />\n'
    )
    f.write(
        '\t\t<text text-anchor="start" x="40" y="150" '
        'style = "font-size: 10; font-weight: bold; fill: #00E408">.25% &#60; Freq &#60; .5%</text>\n\n'
    )
    # key TRIPLE HELICAL
    f.write(
        '\t\t<rect x="295" y="110" width="10" height="10" '
        'style="fill: #00ffe9; stroke-width:1; stroke: #3ca390" />\n'
    )
    f.write(
        '\t\t<text text-anchor="start" x="315" y="120" '
        'style = "font-size: 10; font-weight: bold; fill: #00ffe9">Domains</text>\n\n'
    )
    # set up number line + labels + trip helix
    proteinlength = int(protlength_raw) * 2
    n = 0
    for n in range(0, 4):
        base = 200 + 100 * n
        # triple helical region
        f.write('\n<!-- TRIPLE HELICAL REGION -->\n')
        for domain in domainresults:
            f.write('\n<!-- Domain -->\n')
            f.write(shapes.bluerect(base, domain[0] * 2, domain[1] * 2))

        # line and labels
        f.write('\n<!-- #LINE + LABELS -->\n')
        f.write(
            '\t<line x1="0" y1="' + str(base) + '" x2="' + str(proteinlength) + '" y2="' + str(base) + '" '
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
            '\t<line x1="' + str(proteinlength) + '" y1="' + str(base - 10) + '" x2="' + str(proteinlength) + '" y2="' + str(base + 10) + '" '
                                                                                              'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
        )
        f.write(
            '\t<text text-anchor="end" x="' + str(proteinlength+20) + '" y="' + str(base + 30) + '" '
                                                                        'style = "font-size: 12px;">' + str(proteinlength) + '</text>\n\n'
        )
        # ticks (every 100 AAs)
        f.write('\n<!-- TICKS ON #LINE -->\n')
        position = 0
        while position <= proteinlength:
            f.write(
                '\t<line x1="' + str(position) + '" y1="' + str(base - 5) +
                '" x2="' + str(position) + '" y2="' + str(base + 5) + '" '
                                                                    'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
            )
            position += 200
        n += 1


    # triangle markers
    f.write('\n<!-- MARKERS -->\n')
    for sublist in results:
        num = int(sublist[0])
        fill = ''
        stroke = ''
        if sublist[1] == 2:
            fill = '#FE018C'  # pink
            stroke = '#60002c'
            f.write(shapes.only_triangle(500, num, fill, stroke))
            f.write(shapes.only_triangle(200, num, fill, stroke))

        elif sublist[1] <= 0.0025:
            fill = '#fff600'  # yellow
            stroke = '#593e00'
            f.write(shapes.only_triangle(400, num, fill, stroke))
            f.write(shapes.only_triangle(200, num, fill, stroke))

        elif 0.0025 < sublist[1] <= 0.01:
            fill = '#00E408'  # green
            stroke = '#005908'
            f.write(shapes.only_triangle(300, num, fill, stroke))
            f.write(shapes.only_triangle(200, num, fill, stroke))


    f.write('</svg>\n')