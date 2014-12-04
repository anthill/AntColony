'use strict';

var initRendering = require('./src/rendering.js');
var initializePoints = require('./src/initializePoints.js');
var createEdges = require('./src/createEdges.js');
var initAnts = require('./src/initializeAnts');

module.exports = function init(containerElement, options){

	var render, pointsInfos, edges, population, initVar;


	function _init(containerElement, options){
		pointsInfos = initializePoints(options.nbStart, options.nbRand);
		edges = createEdges(pointsInfos.points);
		population = initAnts(containerElement, pointsInfos, options);
		initVar = {
			pointsInfos: pointsInfos,
			edges: edges,
			population: population
		};
		render = initRendering(containerElement, initVar);
	}

	_init(containerElement, options);

	return {
		togglePlayPause: function(){ render.togglePlayPause() },
		changeOptions: function(opts){
			if (opts.nbStart === undefined && opts.nbRand === undefined){
				// modify population
				render.modifyAnts(opts);
			}
			else{
				// reset previous animation
				render.reset();

				// reset elements
				_init(containerElement, opts);
				
			}
		}
	};
};