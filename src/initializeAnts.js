'use strict'

var antFunction = require('./ant.js');

// var NBANTS = 4000;

module.exports = function (container, options) {

	var Ant = antFunction(container, options);

	var population = new Array(options.nbAnts);
	var possibleStartPointsId = require('./initializePoints.js').possibleStartPointsId;

	for (var i = 0; i < NBANTS; i++) {
	    var newAnt = new Ant(Ant.generateRandStartPoint());
	    newAnt.setDirection();
	    population[i] = newAnt;
	}

	return population;

}