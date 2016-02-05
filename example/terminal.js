// 'use strict'; // This wouldn't allow js to be run in Safari :/

var _antColony = require('../src/index.js');
var inputAnalysis = require('./inputAnalysis.js');
var defaultOptions = require('./defaultOptions.js');

module.exports = function(container){

	var options = defaultOptions;

	// launch animation
	var antColony = _antColony(container, defaultOptions);
	container.addEventListener('click', antColony.togglePlayPause);

	var terminal = document.getElementById('terminal');

	var baseText = '<span class="prompt">ants&gt;&nbsp;</span>';
	var cursor = '<span class="cursor">&nbsp;</span>';
	// var baseTextLength = baseText.length; // 5 is for the <br/>
	var cursorLength = cursor.length;

	var history = [''];
	var currentLine = 0;
	var letterCount = 0;

	function preventDefault(event){

		var key = event.keyCode;

		switch (key){
			case 8: // backspace
				event.preventDefault();
				keyHandler(event);
				break;
			case 13: // enter
				event.preventDefault();
				keyHandler(event);
				break;
			case 32: // spacebar
				event.preventDefault();
				keyHandler(event);
				break;
			case 37: // left
				event.preventDefault();
				break;
			case 38: // up
				event.preventDefault();
				keyHandler(event);
				// terminal.scrollTop -= 19; // 19 should be the actual computed font-size
				break;
			case 39: // right
				event.preventDefault();
				break;
			case 40: // bottom
				event.preventDefault();
				keyHandler(event);
				// terminal.scrollTop += 19; // 19 should be the actual computed font-size
				break;

			default:
				break;
		}
	}

	function keyHandler(event){

		var key = (event.which) ? event.which : 
					((event.charCode) ? event.charCode : 
						((event.keyCode) ? event.keyCode : 0));

		console.log('key pressed ', key);
		// event.preventDefault();

		var line = document.getElementById('currentLine');

		switch (key){
			case 8: // backspace
				if (letterCount > 0) { // check if there are letters to delete
					var text = line.innerHTML.slice(0, - cursorLength);
					line.innerHTML = text.slice(0, -1) + cursor;

					letterCount --;
					// console.log('letter Count :', letterCount);
					// console.log('cursor :', cursor);
					// console.log('baseText :', baseText);
				}
				
				break;
			case 13: // return
				var content = line.innerHTML;
				var length = content.length - cursorLength;

				line.innerHTML = content.slice(0, length);

				// var previousContent = content.slice(0, length);

				var input = content.slice(baseText.length, length); // hit 'enter' seems to add 1 to string
				console.log('input :', input);
				
				history.push(input);
				currentLine = history.length;
				letterCount = 0;
				// console.log('letter Count :', letterCount);
				// console.log('raw :', content.slice(baseTextLength+1, length));
				// console.log('Rcontent :', content + baseText);

				// Check for keywords and launch appropriate function
				var answer = inputAnalysis(input, options, antColony);
				
				var newLine = document.createElement('div');
				newLine.id = 'currentLine';
				newLine.innerHTML = baseText + cursor;
				line.id = '';
			
				if (answer){
					var answerLine = document.createElement('div');
					answerLine.innerHTML = answer;
					terminal.appendChild(answerLine);
				}

				terminal.appendChild(newLine);
				// console.log('newline :', newLine.innerHTML);
				terminal.scrollTop = terminal.scrollHeight;

				break;
			case 38: // up
				var content = line.innerHTML;

				currentLine -= 1;

				if (currentLine < 0)
					currentLine = 0;

				letterCount = history[currentLine].length;
				// console.log('letter Count :', letterCount);

				line.innerHTML = baseText + history[currentLine] + cursor;
				break;

			case 40: // down
				var content = line.innerHTML;

				currentLine += 1;

				// console.log('current ', currentLine);
				// console.log('history length ', history.length);

				if (currentLine >= history.length)
					currentLine = history.length - 1;

				letterCount = history[currentLine].length;
				// console.log('letter Count :', letterCount);

				line.innerHTML = baseText + history[currentLine] + cursor;
				break;
			
			default:
				var content = line.innerHTML;
				var length = content.length - cursorLength;
				// console.log('content :', content);
				var content = content.slice(0, length);
				// console.log('cursor :', cursor);
				// console.log('content :', content);
				letterCount ++;

				line.innerHTML = content + String.fromCharCode(key) + cursor;

				break;
		}
	}

	function activateTerminal(event){
		event.target.focus();
		
		window.addEventListener('keypress', keyHandler);
		window.addEventListener('keydown', preventDefault);
	}

	function deactivateTerminal(event){
		event.target.blur();

		window.removeEventListener('keypress', keyHandler);
		window.removeEventListener('keydown', preventDefault);
	}

	terminal.addEventListener('mouseover', activateTerminal);
	terminal.addEventListener('mouseout', deactivateTerminal);

	return antColony;
};

