function dynamicSize() {
	var windowX = $(window).width();
	var windowY = $(window).height();
	var visX = 900;
	var a_w;
	$("body").css('height', windowY + "px");
	if (windowY > 710) {
		$("#summarybox").css('height', windowY + 'px');
		$('aside').css('height', windowY + 'px');
		$('#container').css('height', windowY + 'px');
		$("#vis").css('height', windowY + 'px');
	} else{
		$("#summarybox").css('height', '710px');
		$('aside').css('height', '710px');
		$('#container').css('height', '710px');
		$("#vis").css('height', '710px');
	}
	if (windowX > 1300) {
		a_w = ((windowX - visX) / 2) - 1.5;
		$("#container").css('width', "100%");
		$("#container").css("height", '100%');
		$("aside").css('width', a_w + 'px');
	} 
	console.log(windowX, windowY);
	if (windowX <= 1300) {
		$("aside").css('width', '198px');
		$("#container").css('width', '1300px');
	}
}

$(window).resize(function() { return dynamicSize(); });

$(window).load(function() { return dynamicSize(); });
