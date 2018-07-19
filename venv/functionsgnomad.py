import re

# parses data
def readresults(dataset):
    resultsraw = []
    pattern = re.compile(
        r'(?<=p\.[A-Z][a-z]{2})\d+(?=[A-Z][a-z]{2}$)'  # searches for p.XXX[Location]XXX"
    )
    intpattern = re.compile(r'[\d.]+(?=e-)')  # finds base in exponent
    intexp = re.compile(r'(?<=e)[-\d]+')  # finds exponent
    glypattern = re.compile(r'(?<=p.Gly)\d+(?=[A-Z][a-z]+)')  # finds exponent
    for row in dataset:
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
            resultsraw.append(result)

    return resultsraw


# checks for ascending
def checkascending(listresults):
    results = []
    errors = []
    minimum = 0
    for result in listresults:
        if result[0] >= minimum:
            minimum = result[0]
            results.append(result)
        elif result[0] < minimum:
            errors.append("Error: %d not ascending from %d" % (result[0], minimum))
    return results, errors


# checks for descending
def checkdescending (listresults):
    results = []
    errors = []
    for index, result in enumerate(listresults):
        if index == 0:
            maximum = result[0]
            results.append(result)
        elif result[0] <= maximum:
            maximum = result[0]
            results.append(result)
        elif result[0] > maximum:
            errors.append("Error: %d not descending from %d" % (result[0], maximum))
    return results, errors

