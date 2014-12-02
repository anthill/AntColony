'use strict';

var initRendering = require('./src/rendering.js');

module.exports = function(containerElement, opts){
	var options = {
		speed: 0.001,
		nbAnts: 4000,
		weight: 10,
		repSize: 0.05,
		repSpeed: 0.002,
		nbStart: 20,
		nbRand: 500
		// obj par defaut
	};

	Object.assign(options, opts);

	if (options.nbStart != 20 || options.nbRand != 500){
		// Load the whole canvas with all options
		initRendering(containerElement, options);
    	var points = require('./src/initializePoints.js')(options.nbStart, options.nbRand);
    	var edges = require('./src/createEdges.js');
	}
	else
	{
		// update directly
		// - speed
		// - nbAnts
		// - weight
		// - repSize
		// - repSpeed

		// get population from 'initAnts'
	}

    
};