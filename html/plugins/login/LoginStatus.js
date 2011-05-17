
dojo.provide("plugins.login.LoginStatus");

// DISPLAY LOGIN STATUS AT RIGHT SIDE OF TOOLBAR

dojo.declare( "plugins.login.LoginStatus",
	[ dijit._Widget, dijit._Templated ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "login/templates/loginstatus.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	constructor: function ()
	{

		this.loadCSS();
	},

	postMixInProperties: function()
	{
	},

	postCreate : function()
	{

		this.startup();
	},

	startup : function ()
	{

		this.inherited(arguments);
	},

	// LOAD CUSTOM STYLE SHEETS
	loadCSS : function ()
	{
		var cssFiles = [ "plugins/login/css/loginstatus.css" ];
		for ( var i in cssFiles )
		{
			var cssFile = cssFiles[i];
			var cssNode = document.createElement('link');
			cssNode.type = 'text/css';
			cssNode.rel = 'stylesheet';
			cssNode.href = cssFile;
			cssNode.media = 'screen';
			//cssNode.title = 'loginCSS';
			document.getElementsByTagName("head")[0].appendChild(cssNode);
		}
	}
}); // end of plugins.login.LoginStatus
