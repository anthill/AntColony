'use strict'

var antFunction = require('./ant.js');

var nbAnts = 3000;

module.exports = function (container) {

	var Ant = antFunction(container);

	var population = new Array(nbAnts);
	var possibleStartPointsId = require('./initializePoints.js').possibleStartPointsId;

	for (var i = 0; i < nbAnts; i++) {
	    var newAnt = new Ant(Ant.generateRandStartPoint());
	    newAnt.setDirection();
	    population[i] = newAnt;
	}

	return population;

}