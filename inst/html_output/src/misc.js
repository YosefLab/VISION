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

function destroyClickedElement(event) {
    document.body.removeChild(event.target);
}

Element.prototype.remove = function() {
    this.parentElement.removeChild(this);
}

