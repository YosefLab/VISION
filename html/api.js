var api = (function(){

    var output = {};

    // Signature API
    
    output.signature = {}

    output.signature.info = function(sig_name){
        var query = "/Signature/Info/"
        query = query.concat(encodeURI(sig_name));
		//return $.ajax(query).then(x => x);
        return $.ajax(query).then(x => JSON.parse(x))
    }

    output.signature.listPrecomputed = function(){
        var query = "/Signature/ListPrecomputed"
		//return $.ajax(query).then(x => x);
        return $.ajax(query).then(x => JSON.parse(x))
    }

    output.signature.scores = function(sig_name){
        var query = "/Signature/Scores/"
        query = query.concat(encodeURI(sig_name))
        //return $.ajax(query).then(x => x);
        return $.ajax(query).then(x => JSON.parse(x))
    }

    output.signature.ranks = function(sig_name){
        var query = "/Signature/Ranks/"
        query = query.concat(encodeURI(sig_name))
        //return $.ajax(query).then(x => x);
        return $.ajax(query).then(x => JSON.parse(x))
    }

    output.signature.expression = function(sig_name){
        var query = "/Signature/Expression/"
        query = query.concat(encodeURI(sig_name))
        //return $.ajax(query).then(x => x);
        return $.ajax(query).then(x => JSON.parse(x))
    }

    output.signature.clusters = function(precomputed) {
		var query = ""
		if (precomputed) {
			query = "/SigClusters/Precomputed"
		} else {
			query = "/SigClusters/Normal"
		}
		//return $.ajax(query).then(x => x);
		return $.ajax(query).then(x => JSON.parse(x))
	}

    // FilterGroup API

    output.filterGroup = {}

    output.filterGroup.sigProjMatrix = function(filter_group, precomputed)
    {
      var query = "/FilterGroup/"
      if (precomputed) {
		query = query.concat(encodeURI(filter_group), "/SigProjMatrix/Precomputed")
	  } else {
	  	query = query.concat(encodeURI(filter_group), "/SigProjMatrix/Normal")
	  }
	  //return $.ajax(query).then(x => x)
      return $.ajax(query).then(x => JSON.parse(x))
    }

    output.filterGroup.sigProjMatrixP = function(filter_group, precomputed)
    {
      var query = "/FilterGroup/"
      if (precomputed) {
		query = query.concat(encodeURI(filter_group), "/SigProjMatrix_P/Precomputed")
	  } else {
	  	query = query.concat(encodeURI(filter_group), "/SigProjMatrix_P/Normal")
	  }
	  //return $.ajax(query).then(x => x);
	  return $.ajax(query).then(x => JSON.parse(x))
    }

    output.filterGroup.list = function()
    {
      var query = "/FilterGroup/list";
      return $.ajax(query).then(x => JSON.parse(x))
    }

    output.filterGroup.genes = function(filter_group)
    {
      var query = "/FilterGroup/"
      query = query.concat(encodeURI(filter_group), "/genes")
      //return $.ajax(query).then(x => x);
      return $.ajax(query).then(x => JSON.parse(x))
    }

    // Projection API

    output.projection = {}
     
    output.projection.coordinates = function(filter_group, projection_name)
    {
      var query = "/FilterGroup/"
      query = query.concat(encodeURI(filter_group), "/",
              encodeURI(projection_name), "/coordinates")
      //return $.ajax(query).then(x => x);
      return $.ajax(query).then(x => JSON.parse(x))
    }

    output.projection.clusters = function(filter_group, projection_name,
            cluster_method, parameter)
    {
      var query = "/FilterGroup/"
      query = query.concat(encodeURI(filter_group), "/", encodeURI(projection_name),
              "/clusters/", encodeURI(cluster_method), "/", encodeURI(parameter))
      //return $.ajax(query).then(x => x);
      return $.ajax(query).then(x => JSON.parse(x))

    }

    output.projection.tree = function(filter_group, projection_name) 
	{
		var query = "/FilterGroup/";
		query = query.concat(encodeURI(filter_group), "/Tree");
		//return $.ajax(query).then(x => x);
		return $.ajax(query).then(x => JSON.parse(x));
	}


    // Expression API


    output.expression = {}

    output.expression.all = function() {
      //return $.ajax(query).then(x => x);
      return $.ajax(query).then(x => JSON.parse(x))
    }


	// Analysis API
	
	output.analysis = {}

	output.analysis.run = function(subset) {
		var query = "/Analysis/Run/";
		return $.ajax({
			type: "POST",
			url: query,
			data: JSON.stringify(subset),
		}).done(alert("Running Subset Analysis"));	
	}

    return output;

})();
