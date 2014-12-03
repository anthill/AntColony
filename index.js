'use strict';

var initRendering = require('./src/rendering.js');

module.exports = function(containerElement, opts){
	var options = {
		velocity: 0.001,
		nbAnts: 4000,
		weight: 10,
		repSize: 0.05,
		repSpeed: 0.002,
		nbStart: 20,
		nbRand: 500
		// obj par defaut
	};

	Object.assign(options, opts);

	if (options.nbStart != 0 || options.nbRand != 0){
		// Load the whole canvas with all options
		initRendering(containerElement, options);
    	var points = require('./src/initializePoints.js')(options.nbStart, options.nbRand);
    	console.log('Nb of points ', points.length);
    	var edges = require('./src/createEdges.js');
    	console.log('Lancé');
	}
	else
	{
		// // Update all ants according to 'options' values
		// var population = require('./src/initializeAnts.js');

		// var deadAnts = -(population.length - options.nbAnts);

		// population = population.slice(deadAnts); // nbAnts

		// population.forEach(function(ant){
		// 	ant.velocity = options.speed;
		// 	ant.weight = options.weight;
		// 	ant.repSize = options.repSize;
		// 	ant.repSpeed = options.repSpeed;
		// });
	}

    
};