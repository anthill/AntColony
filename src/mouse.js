'use strict'

module.exports = function (container){

	var mouse = {
	    x: 0,
	    y: 0
	};

	container.addEventListener( 'mousemove', function(e){
	    var rect = container.getBoundingClientRect();
	    mouse.x = e.clientX / rect.width;
	    mouse.y = e.clientY / rect.height;
	});

	return mouse;

};
