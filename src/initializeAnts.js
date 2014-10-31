'use strict'

var antFunction = require('./ant.js');

var NBANTS = 4000;

module.exports = function (container) {

	var Ant = antFunction(container);

	var population = new Array(NBANTS);
	var possibleStartPointsId = require('./initializePoints.js').possibleStartPointsId;

	for (var i = 0; i < NBANTS; i++) {
	    var newAnt = new Ant(Ant.generateRandStartPoint());
	    newAnt.setDirection();
	    population[i] = newAnt;
	}

	return population;

}