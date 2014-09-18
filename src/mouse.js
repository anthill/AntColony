'use strict'

module.exports = function (container){

	var mouse = {
	    x: 0,
	    y: 0
	};

	container.addEventListener( 'mousemove', function(e){
		console.log(container)
	    var rect = container.getBoundingClientRect();
	    mouse.x = e.clientX / rect.width;
	    mouse.y = e.clientY / rect.height;
	});

	return mouse;

};