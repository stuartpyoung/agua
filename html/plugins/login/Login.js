
dojo.provide("plugins.login.Login");

// DISPLAY A LOGIN DIALOGUE WINDOW AND AUTHENTICATE
// WITH THE SERVER TO RETRIEVE A SESSION ID WHICH IS
// STORED AS Agua.cookie("sessionId")

//dojo.require("plugins.core.Common");

//dojo.require("plugins.admin.Controller");
//dojo.require("plugins.workflow.Controller");
//dojo.require("plugins.project.Controller");
//dojo.require("plugins.report.Controller");
//dojo.require("plugins.help.Controller");



// REQUIRED WIDGETS
dojo.require("dojox.widget.Dialog");
dojo.require("dijit.form.Button");
dojo.require("dijit.form.TextBox");

dojo.require("dojo.fx.easing");
dojo.require("dojox.timing.Sequence");

// HAS A
dojo.require("plugins.login.LoginStatus");

dojo.declare( "plugins.login.Login",
	[ dijit._Widget, dijit._Templated ],
{

	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "login/templates/login.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	loginMessage : "Please enter your CCS password",

	constructor: function ()
	{

		this.loadCSS();

		// GENERATE LOGIN STATUS TABLE IN RIGHT SIDE OF AGUA TOOLBAR
		this.statusBar = new plugins.login.LoginStatus();
		this.statusBar.launcher.innerHTML = "Log In";
		var listener = dojo.connect(this.statusBar.launcher, "onclick", this, "show");
		this.statusBar.launcher.listener = listener;
		Agua.toolbar.domNode.appendChild(this.statusBar.containerNode);
	},

	// postMixInProperties
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

		// ADD CLASS TO DIALOGUE
		dojo.addClass(this.domNode, 'login');

		// SHOW LOGIN WINDOW
		this.show();

		// LINK this.loginDialogue.show TO this.show
		//dojo.connect(this.loginDialogue, "show", this, "show");
	},

	// SHOW 	
	show : function ()
	{

		// RESET STYLE TO DEFAULT
		dojo.removeClass(this.message, "error");
		dojo.removeClass(this.message, "accepted");

		//// GET USERNAME FROM URL, OR USE DEFAULT VALUE
		//var url = window.location.href;
		//var urlUsername = url.match(/\\?(\w+)$/)[0];
		//if ( urlUsername == null || ! urlUsername || urlUsername == "html" )
		//{
		//	urlUsername = '';
		//}
		//this.username.setValue(urlUsername);

		// ADD CLASS
		dojo.addClass(this.loginDialogue.domNode, "login");

		// SET DEFAULT CSS CLASSES
		dojo.removeClass(this.message, "error");
		dojo.removeClass(this.message, "accepted");

		if ( this.loginDialogue == null )
		{
			return;
		}
		this.message.innerHTML = this.loginMessage;

		// FOCUS ON PASSWORD INPUT IF 'RETURN' KEY PRESSED WHILE IN USERNAME INPUT
		var loginObject = this;
		dojo.connect(this.username, "onKeyPress", function(event){
			var key = event.charOrCode;

			// STOP EVENT BUBBLING
			event.stopPropagation();   

			// JUMP TO PASSWORD INPUT IF 'RETURN' KEY PRESSED
			if ( key == 13 )
			{
				loginObject.password.focus();
			}
		});

		// DO LOGIN IF 'RETURN' KEY PRESSED WHILE IN PASSWORD INPUT
		var loginObject = this;
		dojo.connect(this.password, "onKeyPress", function(event){
			var key = event.charOrCode;

			// STOP EVENT BUBBLING
			event.stopPropagation();   

			// LOGIN IF 'RETURN' KEY PRESSED
			if ( key == 13 )
			{
				loginObject.login();
			}

			// QUIT LOGIN WINDOW IF 'ESCAPE' KEY IS PRESSED
			if (key == dojo.keys.ESCAPE)
			{
				// FADE OUT LOGIN WINDOW
				dojo.fadeOut({ node: "loginDialogue", duration: 500 }).play();
				loginObject.loginDialogue.hide();
			}
		});

		this.loginDialogue.show();
	},

	// RESET STYLES TO DEFAULT AND HIDE DIALOGUE
	hide : function ()
	{

		this.loginDialogue.hide();

		dojo.removeClass(this.message, "error");
		dojo.removeClass(this.message, "accepted");

		this.message.innerHTML = this.loginMessage;

	},

	// AUTHENTICATE USER NAME AND PASSWORD
	login : function ()
	{

		var username = this.username.value;
		var password = this.password.attr('value');

		if ( password == null || password == '' )	return;

		dojo.removeClass(this.message, "error");
		dojo.removeClass(this.message, "accepted");
		this.message.innerHTML = "Authenticating...";

		// SET SESSION ID AND USER NAME TO NULL	
		Agua.cookie("sessionId", null); 
		Agua.cookie("username", null); 

		// CREATE JSON QUERY
		var query = new Object;
		query.username = username;
		query.password = password;
		query.mode = "login";

		var url = Agua.cgiUrl + "agua?";


		var loginObject = this;


		// DO xhrPut TO SERVER
		var deferred = dojo.xhrPut(
			{
				url: url,
				contentType: "json",
				handleAs: 'json',
				sync: true,
				putData: dojo.toJson(query),
				preventCache: true,
				timeout: 20000,
				load: function(data) {


					// SHOW ERROR IF PRESENT AS data.error
					if ( data.error != null )
					{

						// SHOW ERROR STATUS ON LOGIN BUTTON
						dojo.addClass(loginObject.message, "error");
						loginObject.message.innerHTML = data.error;
					}

					// IF NO ERROR, PROCESS LOGIN
					else
					{

						// SET sessionId AND username IN DOJO COOKIE
						// THE COOKIE WILL EXPIRE AT THE END OF THE SESSION
						Agua.cookie("sessionId", data.sessionId);
						Agua.cookie("username", username);


						// SHOW 'Accepted' IN GREEN
						dojo.removeClass(loginObject.message, "error");
						dojo.addClass(loginObject.message, "accepted");
						loginObject.message.innerHTML = "Accepted";

						// CHANGE LOGIN STATUS BAR TO 'username' AND 'Log Out'
						loginObject.statusBar.username.innerHTML = username;
						loginObject.statusBar.launcher.innerHTML = "Log Out";

						// REMOVE EXISTING 'login' LISTENER AND CREATE NEW 'logout' LISTENER
						dojo.disconnect(loginObject.statusBar.launcher.listener);
						var listener = dojo.connect(loginObject.statusBar.launcher, "onclick", loginObject, "logout");
						loginObject.statusBar.launcher.listener = listener;


						//// FADE OUT LOGIN WINDOW
						//loginObject.hide();
						//

						////// INITIALISE PLUGINS
						//setTimeout( function(thisObj){ Agua.startPlugins(); }, 100, loginObject);

						// FADE OUT LOGIN WINDOW
						//dojo.fadeOut({ node: loginObject.loginDialogue, duration: 500 }).play();
						setTimeout( function(thisObj){ thisObj.hide(); }, 100, loginObject);
					}					
				},
				error: function(response, ioArgs) {
					// NO ERROR MESSAGE TO AVOID THIS:
					// Login.Error with JSON Post, response: {"message":"this.loginButton is undefined"
					return response;
				}
			}
		);

	},


	handleLogin : function (response)
	{

	},


	logout : function (e)
	{
		// SET COOKIE sessionId TO NULL
		Agua.cookie('sessionId', null);
		Agua.cookie('username', null);

		// RESET LOGIN DIALOGUE
		this.message.innerHTML = this.loginMessage;
		dojo.removeClass(this.message, "error");
		dojo.removeClass(this.message, "accepted");

		// RESET LOGIN STATUS
		this.statusBar.username.innerHTML = 'Logged out';
		this.statusBar.launcher.innerHTML = 'Log In';

		// CREATE NEW LISTENER
		dojo.disconnect(this.statusBar.launcher.listener);
		var listener = dojo.connect(this.statusBar.launcher, "onclick", this, "show");
		this.statusBar.launcher.listener = listener;
	},


	// LOAD CUSTOM STYLE SHEETS
	loadCSS : function ()
	{
		var cssFiles = [ "plugins/login/css/login.css" , "dojo-1.5.0/dojox/widget/Dialog/Dialog.css" ];
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
}); // end of plugins.login.Login
