# The main readme.md at the repository's root is used as the landing
# page of the doxygen documentation, which becomes an issue when
# resolving links (github vs local docs). This script creates a
# separate readme.md for doxygen and replaces the links to point to
# local or github links, depending in the provided arguments.

# --- STD Imports ---
import pathlib
import re
import argparse


parser = argparse.ArgumentParser("prepareDocs")
parser.add_argument("-w",
                    dest = "web",
                    action = "store_const",
                    const = True,
                    default = False,
                    help = "prepare docs for the web")
arguments = parser.parse_args()

rootPath = "https://matekelemen.github.io/mcgs/" if arguments.web else "../../"
docPath = rootPath + "docs/html/"
assetPath = rootPath + ".github/assets/"

docLinkPattern = re.compile(r"\[mcgs::(\w+)\]\(([\w0-9\.#]+)\)")
docLinkReplace = r"[mcgs::\1](" + docPath + r"\2)"

assetLinkPattern = re.compile(r"\.github/assets/")
assetLinkReplace = assetPath

formulaPattern = re.compile(r"\$(.*)\$")
formulaReplace = r"@f$\1@f$"

readme = ""
with open(pathlib.Path(__file__).absolute().parent.parent / "readme.md", "r") as inputFile:
    readme = inputFile.read()

for pattern, replace in ((docLinkPattern, docLinkReplace), (assetLinkPattern, assetLinkReplace), (formulaPattern, formulaReplace)):
    readme = re.sub(pattern, replace, readme)

with open(pathlib.Path(__file__).parent / "readme.md", "w") as outputFile:
    outputFile.write(readme)
