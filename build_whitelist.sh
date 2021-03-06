# This file is used to build a whitelist
# of files that the FastProject Viewer is
# allowed to serve.

git ls-files | \
    grep "inst/html_output" | \
    egrep "(.js|.html|.css|.png|.jpg|.gif|.svg)$" | \
    sed "s.inst/html_output/.html_output/." > \
    inst/html_output/whitelist.txt
