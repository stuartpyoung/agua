dojo.provide("plugins.login.Controller");

// OBJECT:  plugins.login.Controller

var loginController;

// HAS
dojo.require("plugins.login.Login");
dojo.require("plugins.login.LoginStatus");

dojo.declare( "plugins.login.Controller", 
	[],
{
	// PANE ID COUNTER (INCREMENTS FOR EACH NEW WORKFLOW TAB CREATED)
	paneId : 1,

	login : null,


	// CONSTRUCTOR	
	constructor : function(args) {


		this.loadCSS();
	},


	loadCSS : function()
	{
		// LOAD CSS
		var cssFiles = [ "plugins/login/css/login.css" ];
		for ( var i in cssFiles )
		{
			var cssFile = cssFiles[i];
			var cssNode = document.createElement('link');
			cssNode.type = 'text/css';
			cssNode.rel = 'stylesheet';
			cssNode.href = cssFile;
			cssNode.media = 'screen';
			document.getElementsByTagName("head")[0].appendChild(cssNode);
		}
	}

}); // end of Controller

dojo.addOnLoad(
	function()
	{

		//Agua.loginController = new plugins.login.Login();

		// COMMENT THIS TO DISABLE AUTOMATIC LOGIN
		//Agua.loginController.login();
	}
);
