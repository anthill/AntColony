'use strict';

var _antColony = require('../index.js');

var container = document.querySelector('.colony');

var options = {
	velocity: 0.001,
	nbAnts: 3000,
	intelligence: 0.95,
	repSize: 0.05,
	repSpeed: 0.002,
	nbStart: 500,
	nbRand: 1000
	// obj par defaut
};

var antColony = _antColony(container, options);

window.addEventListener('click', function (){
	// options.velocity = 0.003;
	options.nbAnts = 3000;
	// options.weight = 10000000;
	// options.repSpeed = 0.01;
	// options.repSize = 0.1;

	// antColony.changeOptions(options);
	antColony.changeOptions(options);
});

