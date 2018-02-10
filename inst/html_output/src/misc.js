function createSigModal(sig_key){

    return api.signature.info(sig_key).then(function(sig_info){

        var sig_data = [];
        for(var gene in sig_info['sigDict'])
        {
            sig_data.push({'Gene': gene, 'Sign': sig_info['sigDict'][gene]});
        }
        var sigModal = $('#signatureModal');
        sigModal.find('h4').text(sig_key);
        var tableRows = d3.select('#signatureModal').select('tbody').selectAll('tr')
            .data(sig_data);
        tableRows.enter().append('tr');
        tableRows.exit().remove();

        var tableCells = tableRows.selectAll('td').data(function(d){return [d.Gene, d.Sign];});
        tableCells.enter().append('td');
        tableCells.text(function(d, i){
            if(i == 0){return d;}
            else{
                if(d == 1)  {return "+";}
                if(d == -1) {return "-";}
                if(d == 0)  {return "Unsigned";}
                return "Unknown";
            }
        });

        tableCells.exit().remove();

        sigModal.modal();
    });
}



function createGeneModal()
{
    return api.filterGroup.genes(global_status.filter_group)
        .then(function(genes){

            //Calculate max width
            var width_and_index = genes.map(function(e,i){return [e.length, i]});
            width_and_index.sort(function(a,b){return Math.sign(b[0] - a[0]);});
            var top10 = width_and_index.slice(0,10).map(function(e){return genes[e[1]];});
            var widths = [];
            for(var i = 0; i < top10.length; i++)
            {
                var div = document.createElement("div");
                $(div).text(top10[i]).css("position","absolute").css("left", "-9999px");
                $('body').append(div);
                widths.push($(div).width());
            }

            var maxWidth = d3.max(widths);

            /* var geneDivs = d3.select('#geneModal').select('.modal-body').selectAll('di
                .data(genes.sort());
                */

            geneDivs.enter().append('div');
            geneDivs.exit().remove();

            geneDivs
                .text(function(d){return d;})
                .style("width", maxWidth + "px");

            $('#geneModal').modal();
        });
}


function exprotSigProj() {

    var sig_key = global_status.plotted_signature;
    var proj_key = global_status.plotted_projection;
    var filter_group = global_status.filter_group;

    if (sig_key.length == 0 && proj_key.length == 0) {
        $("#plot_title_div").children().eq(0).text("");
        $("#plot_title_div").children().eq(1).text("");
        global_scatter.setData([], false);
        return $().promise();
    }

    var proj_promise = api.projection.coordinates(filter_group, proj_key);

    var sig_promise;
    if (global_status.scatterColorOption == "value" || global_data.sigIsMeta[sig_key]) {
        sig_promise = api.signature.scores(sig_key)
    } else if (global_status.scatterColorOption == "rank") {
        sig_promise = api.signature.ranks(sig_key)
    }

    return $.when(proj_promise, sig_promise)
        .then(function(projection, signature) {

            var points = [];
            for (var sample_label in signature) {
                var x = projection[sample_label][0];
                var y = projection[sample_label][1];
                var sig_score = signature[sample_label];
                points.push([x, y, sig_score, sample_label]);
            } 

            var lineArray = [];
            points.forEach(function(infoArray, index) {
                var line = inforArray.join(",");
                lineArray.push(index == 0 ? "data:text/csv;charset=utf-8," + line : line);
            });
            var csvContent = lineArray.join("\n");

            var encodeURI = encodeURI(csvContent);
            var downloadLink = document.createElement("a");
            downloadLink.setAttribute("href", encodedURI);
            downloadLink.setAttribute("download", proj_key + ".csv");
            downloadLink.onclick = destroyClickedElement;
            downloadLink.style.display = "none";
            document.body.appendChild(downloadLink);

            downloadLink.click();
        });

}


function destroyClickedElement(event) {
    document.body.removeChild(event.target);
}

Element.prototype.remove = function() {
    this.parentElement.removeChild(this);
}

