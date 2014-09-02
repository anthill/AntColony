'use strict'

var Ant = require('./ant.js');

var nbAnts = 500;

var population = new Array(nbAnts);
var possibleStartPointsId = require('./initializePoints.js').possibleStartPointsId;

for (var i = 0; i < nbAnts; i++) {
    var newAnt = new Ant(Ant.generateRandStartPoint());
    newAnt.setDirection();
    population[i] = newAnt;
}

module.exports = population;