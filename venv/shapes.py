def bluerect(height, start, end):
    return '\t<polygon points="' + str(start) + ',' + str(height - 6) + ' ' \
           + str(end) + ',' + str(height - 6) + ' ' + str(end) + ',' \
           + str(height + 6) + ' ' + str(start) + ',' + str(height + 6) \
           + '" style = "fill: #9bfffb; fill-opacity: 0.6; stroke: #22e2db; stroke-width:1" />\n\n'


def greenrect( height, start, end):
    return '\t<polygon points="' + str(start) + ',' + str(height - 6) + ' ' + str(end) + ',' + str(height - 6) + ' ' \
           + str(end) + ',' + str(height + 6) + ' ' + str(start) + ',' + str(height + 6) + '" ' \
            'style = "fill:#b6ff9b; fill-opacity: 0.6; stroke: #00E408; stroke-width:1" />\n\n'


def triangle(height, position, fill, stroke):
    start = str(position * 3 + 1 - 7)
    mid = str(position * 3 + 1)
    end = str(position * 3 + 1 + 7)
    return  '\t<polygon points="' + start + ','+ str(height-30) + ' ' + end + ',' + str(height-30) + ' ' + mid + ',' \
            + str(height-15) + '" style = "fill: ' + fill + '; stroke: ' + stroke + '; stroke-width:1" />\n'


def line(height, position, stroke):
    mid = str(position * 3 + 1)
    return '\t<line x1="' + mid + '" y1="' + str(height - 16) + '" x2="' + mid + '" y2="' + str(height - 1) + '" ' \
    'style="stroke: ' + stroke + '; stroke-width:2" />\n'


def circle(height, position):
    return '<circle cx="' + position + '" cy="' + str(height - 38) + '" r="4" stroke-width=".5" /> \n'