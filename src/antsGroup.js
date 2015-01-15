'use strict'

// var antFunction = require('./ant.js');

// var NBANTS = 4000;

module.exports = function (Ant) {

	// var Ant = antFunction(container, pointsInfos, options);
	var nbAntsPerStep = 100;

	// var population = new Array(options.nbAnts);
	// var possibleStartPointsId = pointsInfos.possibleStartPointsId;

	function createGroup(population){
		for (var i = 0; i < nbAntsPerStep; i++) {
			var newAnt = new Ant(Ant.generateRandStartPoint());
			newAnt.setDirection();
			population.push(newAnt);
		}

		console.log('Created Ants Group: \
(+ ' + nbAntsPerStep + ') => ' + population.length);

		return population;
	}

	function removeGroup(population, nbDead){
		population = population.slice(0, population.length - nbDead);

		console.log('Removed Ants Group: \
(- ' + nbAntsPerStep + ') => ' + population.length);

		return population;

	}

	return {
		create: createGroup,
		remove: removeGroup
	};

}
	