var api = (function(){

    var output = {};

    // Signature API
    
    output.signature = {}

    output.signature.info = function(sig_name){
        var query = "/Signature/Info/"
        query = query.concat(encodeURI(sig_name))
        return $.ajax(query)
    }

    output.signature.scores = function(sig_name){
        var query = "/Signature/Scores/"
        query = query.concat(encodeURI(sig_name))
        return $.ajax(query)
    }

    output.signature.ranks = function(sig_name){
        var query = "/Signature/Ranks/"
        query = query.concat(encodeURI(sig_name))
        return $.ajax(query)
    }

    // FilterGroup API

    output.filterGroup = {}

    output.filterGroup.sigProjMatrix = function(filter_group)
    {
      var query = "/FilterGroup/"
      query = query.concat(encodeURI(filter_group), "/SigProjMatrix")
      return $.ajax(query)
    }

    output.filterGroup.sigProjMatrixP = function(filter_group)
    {
      var query = "/FilterGroup/"
      query = query.concat(encodeURI(filter_group), "/SigProjMatrix_P")
      return $.ajax(query)
    }

    // Projection API

    output.projection = {}
     
    output.projection.coordinates = function(filter_group, projection_name)
    {
      var query = "/FilterGroup/"
      query = query.concat(encodeURI(filter_group), "/",
              encodeURI(projection_name), "/coordinates")
      return $.ajax(query)
    }

    output.projection.clusters = function(filter_group, projection_name,
            cluster_method)
    {
      var query = "/FilterGroup/"
      query = query.concat(encodeURI(filter_group), "/", encodeURI(projection_name),
              "/clusters/", encodeURI(cluster_method))
      return $.ajax(query)

    }


    // Expression API


    output.expression = {}

    output.expression.all = function() {
      return $.ajax('/Expression')
    }

    return output;

})();
