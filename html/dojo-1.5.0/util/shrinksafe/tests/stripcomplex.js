result = "";

/* For the purpose of these tests, we have no actual console with
   which to test output, so create a dummy console and capture to
   'result', which is examined by the test framework post eval. */
var console = {
	log: function(arg) { result += "LOG: " + arg; },
	debug: function(arg) { result += "DEBUG: " + arg; },
	warn: function(arg) { result += "WARN: " + arg; },
	error: function(arg) { result += "ERROR: " + arg; },
	dir: function(arg) { result += "DIR: " + arg; }
};
(function() {
	var fn = function(arg) {
		return "fn saw arg '" + arg + "'.";
	}
	var str = "wooo";
	console.debug(str + "some \\ dodgy \" characters *$!?//" + (1+2) + (~2) + fn(str));
	console.error(str + "some \\ dodgy \" characters *$!?//" + (1+2) + (~2) + fn(str));

	// from ticket http://bugs.dojotoolkit.org/ticket/8549

	//It would be interesting to see how this comes out:
	if(true)
	else
	   var two = "two";

	var bar;

	(function() {
		return console.debug("Debug return statement.");
	})();
	(function() {
		return console.error("Error return statement.");
	})();
})();