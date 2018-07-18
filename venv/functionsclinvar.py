import re


def readresultsascending(table, gene_name, protein_name):
    resultsraw = []
    just_location = []
    errorsraw = []
    location_pat = re.compile(
        r'(?<=NP_\d{6}\.\d:p\.[A-Z][a-z][a-z])\d+(?=[A-Z][a-z]{2}$)'
        # searches for NP_######.#:p.XXX[Location]XXX"
    )
    protein_pat = re.compile(r'(?<=NP_)\d{6}')  # searches for NP_[protein #]
    glypattern = re.compile(r'(?<=p.Gly)\d+(?=[A-Z][a-z]+)')
    for row in table:
        gene = row["symbol"]
        pathogenic = row["pathogenic"]
        likely_path = row["likely_pathogenic"]
        uncertain = row["uncertain_significance"]
        consequence = row["hgvs_p"]

        if gene == gene_name and int(uncertain) >= 1 and location_pat.findall(consequence):
            protein = str(protein_pat.findall(consequence))[2:-2]
            if protein == protein_name:
                key = 0
                gly = 0
                location = str(location_pat.findall(consequence))[2:-2]
                location = int(location)
                if glypattern.findall(consequence):
                    gly = 1
                result = [location, key, gly]
                resultsraw.append(result)
            else:
                errorsraw.append("Error: protein " + protein + " not protein 001839")

        elif gene == gene_name and int(likely_path) >= 1 and location_pat.findall(consequence):
            protein = str(protein_pat.findall(consequence))[2:-2]
            if protein == protein_name:
                location = str(location_pat.findall(consequence))[2:-2]
                location = int(location)
                key = 1
                gly = 0
                if glypattern.findall(consequence):
                    gly = 1
                result = [location, key, gly]
                resultsraw.append(result)
            else:
                errorsraw.append("Error: protein " + protein + " not protein 001839")

        elif gene == gene_name and int(pathogenic) >= 1 and location_pat.findall(consequence):
            protein = str(protein_pat.findall(consequence))[2:-2]
            if protein == protein_name:
                location = str(location_pat.findall(consequence))[2:-2]
                location = int(location)
                key = 2
                gly = 0
                if glypattern.findall(consequence):
                    gly = 1
                result = [location, key, gly]
                resultsraw.append(result)
            else:
                errorsraw.append("Error: protein " + protein + " not protein 001839")
    for result in resultsraw:
        just_location.append(result[0])  # appends first # in each pair to just_location
    for index, item in enumerate(just_location):  # checks if in ascending order, if not print error msg
        if item < just_location[index - 1] and index != 0:
            errorsraw.append("Error at row %d: %d not ascending" % (index, item))
    return errorsraw, resultsraw

def readresultsdescending(table, gene_name, protein_name):
    resultsraw = []
    just_location = []
    errorsraw = []
    location_pat = re.compile(
        r'(?<=NP_\d{6}\.\d:p\.[A-Z][a-z][a-z])\d+(?=[A-Z][a-z]{2}$)'
        # searches for NP_######.#:p.XXX[Location]XXX"
    )
    protein_pat = re.compile(r'(?<=NP_)\d{6}')  # searches for NP_[protein #]
    glypattern = re.compile(r'(?<=p.Gly)\d+(?=[A-Z][a-z]+)')
    for row in table:
        gene = row["symbol"]
        pathogenic = row["pathogenic"]
        likely_path = row["likely_pathogenic"]
        uncertain = row["uncertain_significance"]
        consequence = row["hgvs_p"]

        if gene == gene_name and int(uncertain) >= 1 and location_pat.findall(consequence):
            protein = str(protein_pat.findall(consequence))[2:-2]
            if protein == protein_name:
                key = 0
                gly = 0
                location = str(location_pat.findall(consequence))[2:-2]
                location = int(location)
                if glypattern.findall(consequence):
                    gly = 1
                result = [location, key, gly]
                resultsraw.append(result)
            else:
                errorsraw.append("Error: protein " + protein + " not protein 001839")

        elif gene == gene_name and int(likely_path) >= 1 and location_pat.findall(consequence):
            protein = str(protein_pat.findall(consequence))[2:-2]
            if protein == protein_name:
                location = str(location_pat.findall(consequence))[2:-2]
                location = int(location)
                key = 1
                gly = 0
                if glypattern.findall(consequence):
                    gly = 1
                result = [location, key, gly]
                resultsraw.append(result)
            else:
                errorsraw.append("Error: protein " + protein + " not protein 001839")

        elif gene == gene_name and int(pathogenic) >= 1 and location_pat.findall(consequence):
            protein = str(protein_pat.findall(consequence))[2:-2]
            if protein == protein_name:
                location = str(location_pat.findall(consequence))[2:-2]
                location = int(location)
                key = 2
                gly = 0
                if glypattern.findall(consequence):
                    gly = 1
                result = [location, key, gly]
                resultsraw.append(result)
            else:
                errorsraw.append("Error: protein " + protein + " not protein 001839")
    for result in resultsraw:
        just_location.append(result[0])  # appends first # in each pair to just_location
    for index, item in enumerate(just_location):  # checks if in ascending order, if not print error msg
        if item > just_location[index - 1] and index != 0:
            errorsraw.append("Error at row %d: %d not descending" % (index, item))
    return errorsraw, resultsraw


def readresultsfkrp(table, gene_name, protein_name):
    resultsraw = []
    just_location = []
    errorsraw = []
    location_pat = re.compile(
        r'(?<=NP_\d{6}\.\d:p\.[A-Z][a-z][a-z])\d+(?=[A-Z][a-z]{2}$)'
        # searches for NP_######.#:p.XXX[Location]XXX"
    )
    protein_pat = re.compile(r'(?<=NP_)\d{6}')  # searches for NP_[protein #]
    terpattern = re.compile(r'(?<=p.[A-Z][a-z]{2})\d+(?=Ter)')
    for row in table:
        gene = row["symbol"]
        pathogenic = row["pathogenic"]
        likely_path = row["likely_pathogenic"]
        uncertain = row["uncertain_significance"]
        consequence = row["hgvs_p"]

        if gene == gene_name and int(uncertain) >= 1 and location_pat.findall(consequence):
            protein = str(protein_pat.findall(consequence))[2:-2]
            if protein == protein_name:
                key = 0
                location = str(location_pat.findall(consequence))[2:-2]
                location = int(location)
                if terpattern.findall(consequence):
                    ter = 1
                result = [location, key, ter]
                resultsraw.append(result)
            else:
                errorsraw.append("Error: protein " + protein + " not protein 001839")

        elif gene == gene_name and int(likely_path) >= 1 and location_pat.findall(consequence):
            protein = str(protein_pat.findall(consequence))[2:-2]
            if protein == protein_name:
                location = str(location_pat.findall(consequence))[2:-2]
                location = int(location)
                key = 1
                if terpattern.findall(consequence):
                    ter = 1
                result = [location, key, ter]
                resultsraw.append(result)
            else:
                errorsraw.append("Error: protein " + protein + " not protein 001839")

        elif gene == gene_name and int(pathogenic) >= 1 and location_pat.findall(consequence):
            protein = str(protein_pat.findall(consequence))[2:-2]
            if protein == protein_name:
                location = str(location_pat.findall(consequence))[2:-2]
                location = int(location)
                key = 2
                if terpattern.findall(consequence):
                    ter = 1
                result = [location, key, ter]
                resultsraw.append(result)
            else:
                errorsraw.append("Error: protein " + protein + " not protein 001839")
    for result in resultsraw:
        just_location.append(result[0])  # appends first # in each pair to just_location
    for index, item in enumerate(just_location):  # checks if in ascending order, if not print error msg
        if item < just_location[index - 1] and index != 0:
            errorsraw.append("Error at row %d: %d not ascending" % (index, item))
    return errorsraw, resultsraw


def percent_pathogenic(total, helix_start, helix_end):
    numerator = 0
    denom = 0
    for value in total:
        if helix_start <= value[0] <= helix_end and (value[1] == 1 or value[1] == 2):
            numerator += 1
        if value[1] == 1 or value[1] == 2:
            denom += 1
    percentage = (float(numerator) / float(denom)) * 100
    return numerator, denom, percentage


def percent_vus(total, helix_start, helix_end):
    numerator1 = 0
    denom1 = 0
    for value in total:
        if helix_start <= value[0] <= helix_end and value[1] == 0:
            numerator1 += 1
        if value[1] == 0:
            denom1 += 1
    percentage1 = (float(numerator1) / float(denom1)) * 100
    return numerator1, denom1, percentage1

