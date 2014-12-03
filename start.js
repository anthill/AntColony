'use strict';

var start = require('./index.js');

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

window.addEventListener('click', function (){
	options.velocity = 0.005;
	start(container, options)
});


start(container, options);

