'use strict'

var antFunction = require('./ant.js');

// var NBANTS = 4000;

module.exports = function (container, pointsInfos, options) {

	var Ant = antFunction(container, pointsInfos, options);

	var population = new Array(options.nbAnts);
	var possibleStartPointsId = pointsInfos.possibleStartPointsId;

	for (var i = 0; i < options.nbAnts; i++) {
	    var newAnt = new Ant(Ant.generateRandStartPoint());
	    newAnt.setDirection();
	    population[i] = newAnt;
	}

	return population;

}