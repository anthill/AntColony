'use strict';

var _antColony = require('./index.js');

var container = document.querySelector('.colony');

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

var antColony = _antColony(container, options);

window.addEventListener('click', function (){
	options.nbStart = 0;
	options.nbRand = 0;

	// options.velocity = 0.003;
	// options.nbAnts = 1000;
	// options.weight = 10000000;
	options.repSpeed = 0.01;
	options.repSize = 0.1;

	antColony.changeOptions(options);
});

