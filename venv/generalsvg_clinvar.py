import re
import csv
import os
import argparse
import functionsclinvar as fc
import shapes

parser = argparse.ArgumentParser(description="sets file name")
parser.add_argument("startfile", help="Enter name of data file")
parser.add_argument("destination", help="Enter name of destination file")
args = parser.parse_args()

if os.path.exists(args.startfile) and os.path.exists(args.destination):
    print "All files exist--running program"

    # completes url for accessing pfam/uniprot
    proteinid = raw_input('Enter Protein ID: ')
    gene_name = fc.genename(proteinid)

    # isoform
    accessionid = raw_input('Enter Isoform Accession # (NP_###): ')

    # for checking ascending/descending later
    order = raw_input('Ascending or Descending? (enter A or D): ')
    while order != 'A' and order != 'D' and order != 'a' and order != 'd':
        print "Error: please enter A or D"
        order = raw_input('Ascending or Descending? (enter A or D): ')

    # domains
    repeats = raw_input('How many kinds of domains?: ')
    while repeats.isdigit() is False:
        print "Error: not a number"
        repeats = raw_input('How many kinds of domains? ')
    # prints out info for check
    domainresults, lengthprotein = fc.domains(repeats, proteinid)
    print 'Gene: ' + gene_name
    print 'Length: ' + str(lengthprotein) + ' amino acids'

    # editable bases
    edits_list = []
    edits = raw_input("Would you like to include editable bases? (enter Y or N): ")
    while edits != "Y" and edits != "N" and edits != "y" and edits != "n":
        print "Error: please enter Y or N"
        edits = raw_input('Would you like to include editable bases? (enter Y or N): ')
    if edits == "Y" or edits == "y":
        editfile = raw_input("Enter file directory: ")
        if os.path.exists(editfile):
            with open(editfile) as tsv:
                editsfile = csv.reader(tsv, delimiter="\t")
                edits_list = fc.edits(editsfile, accessionid, len(str(accessionid)))
        while os.path.exists(editfile) is False:
            print "Error: please enter valid directory"
            editfile = raw_input("Enter file directory: ")
    print("Gathering data...")

    # parse data
    results = []
    with open(args.startfile) as tsv:
        reader = csv.DictReader(tsv, delimiter="\t")
        resultsraw, errors = fc.readresults(reader, gene_name, accessionid)
        if order == 'A' or order == 'a':
            results, errors = fc.checkascending(resultsraw, errors)
        elif order == 'D' or order == 'd':
            results, errors = fc.checkdescending(resultsraw, errors)
        if len(results) == 0:
            print "Error: no results found"
            raise SystemExit

        else:

            print 'Writing file...'
            with open(args.destination, 'w') as f:

                for error in errors:
                    f.write('<!--' + error + '-->\n')

                f.write('<svg height="500" width="' + str(lengthprotein * 3 + 20)
                        + '" xmlns="http://www.w3.org/2000/svg">\n\n')

                n = 0
                for n in range(0, 2):
                    height = 50 + 100 * n
                    length = lengthprotein * 3

                    # Domain
                    for domain in domainresults:
                        f.write('\n<!-- Domain -->\n')
                        f.write(shapes.bluerect(height, domain[0]*3+1, domain[1]*3+1))

                    # line and labels
                    f.write('\n<!-- #LINE + LABELS -->\n')
                    f.write(
                        '\t<line x1="0" y1="' + str(height) + '" x2="' + str(length) + '" y2="' + str(height) + '" '
                        'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
                    )
                    f.write(
                        '\t<line x1="1" y1="' + str(height - 10) + '" x2="1" y2="' + str(height + 10) + '" '
                        'style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
                    )
                    f.write(
                        '<text text-anchor="start" x="0" y="' + str(height + 30) + '" '
                        'style = "font-size: 12px;">0</text>'
                    )
                    f.write(
                        '\t<line x1="' + str(length + 2) + '" y1="' + str(height - 10)
                        + '" x2="' + str(length+2) + '" y2="' + str(height + 10)
                        + '" style="stroke: #000000; stroke-width:2; stroke-linecap: round" />\n'
                    )
                    f.write(
                        '\t<text text-anchor="end" x="' + str(length + 10) + '" y="' + str(height + 30) + '" '
                        'style = "font-size: 12px;">' + str(length / 3) + '</text>\n\n'
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
                    n += 1

                # triangle markers
                f.write('\n<!-- MARKERS -->\n')
                for sublist in results:
                    num = sublist[0]
                    fill = ""
                    stroke = ""

                    # pathogenic
                    if sublist[1] == 2 or sublist[1] == 1:
                        stroke = '#ff5e6b'  # pink
                        fill = '#ffa3aa'
                        f.write(shapes.triangle(150, num, fill, stroke))
                        f.write(shapes.line(150, num, stroke))

                    # uncertain
                    elif sublist[1] == 0:
                        fill = '#c2ffa3'  # green
                        stroke = '#478c24'

                    # nonsense
                    if (sublist[1] == 1 or sublist[1] == 2) and sublist[2] == 2:
                        fill = '#ffde72'  # yellow
                        stroke = '#ffa716'
                        f.write(shapes.triangle(150, num, fill, stroke))
                        f.write(shapes.line(150, num, stroke))

                    # gly
                    elif (sublist[1] == 1 or sublist[1] == 2) and sublist[2] == 1:
                        fill = '#9f9eff'  # purple
                        stroke = '#7270ff'
                        f.write(shapes.triangle(150, num, fill, stroke))
                        f.write(shapes.line(150, num, stroke))

                    f.write(shapes.triangle(50, num, fill, stroke))
                    f.write(shapes.line(50, num, stroke))

                # editable bases
                for sublist in edits_list:
                    mid1 = str(int(sublist) * 3 + 1)
                    f.write(shapes.circle(50, mid1))
                    f.write(shapes.circle(150, mid1))

                f.write('</svg>')
else:
    print "error"


# VUS strip
# path
# editable
# founder mutations for gne
# dysf -- mutations we have fibroblasts for
