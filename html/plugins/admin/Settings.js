dojo.provide("plugins.admin.Settings");

// ALLOW THE ADMIN USER TO ADD, REMOVE AND MODIFY USERS

// NEW USERS MUST HAVE username AND email

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dijit.form.Button");
dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.Textarea");
dojo.require("dojo.parser");
dojo.require("dojo.dnd.Source");
dojo.require("dijit.form.ValidationTextBox");
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.admin.UserRow");



dojo.declare(
    "plugins.admin.ValidationTextBox",
	[ dijit.form.ValidationTextBox ],
{

validate: function(/*Boolean*/ isFocused){
	// summary:
	//		Called by oninit, onblur, and onkeypress.
	// description:
	//		Show missing or invalid messages if appropriate, and highlight textbox field.
	// tags:
	//		protected

	// SKIP VALIDATE WHEN LOADING WIDGET
	if ( this.parentWidget == null )	return;

	// IF this IS THE newPassword WIDGET, RUN VALIDATE ON
	// ITS TARGET: THE confirmPassword WIDGET
	if ( this.target != null )
	{
		return this.target.validate(isFocused);
	}
	var message = "";
	var isValid = this.parentWidget.passwordsMatch();

	if(isValid){ this._maskValidSubsetError = true; }
	var isEmpty = this._isEmpty(this.textbox.value);
	var isValidSubset = !isValid && !isEmpty && isFocused && this._isValidSubset();
	this.state = ((isValid || ((!this._hasBeenBlurred || isFocused) && isEmpty) || isValidSubset) && this._maskValidSubsetError) ? "" : "Error";
	if(this.state == "Error"){ this._maskValidSubsetError = isFocused; } // we want the error to show up afer a blur and refocus

	this._setStateClass();
	dijit.setWaiState(this.focusNode, "invalid", isValid ? "false" : "true");
	if(isFocused){
		if(this.state == "Error"){
			message = this.getErrorMessage(true);
		}else{
			message = this.getPromptMessage(true); // show the prompt whever there's no error
		}
		this._maskValidSubsetError = true; // since we're focused, always mask warnings
	}

	this.displayMessage(message);


	return isValid;
}	


});




dojo.declare(
    "plugins.admin.Settings",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "admin/templates/settings.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

//addingUser STATE
addingUser : false,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ "plugins/admin/css/settings.css" ],

// PARENT WIDGET
parentWidget : null,


/////}}} 
constructor : function(args)  {
	// GET INFO FROM ARGS
	this.parentWidget = args.parentWidget;
	this.settings = args.parentWidget.settings;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();		
},

postCreate : function() {

	this.startup();
},


startup : function () {

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 

	// ADD ADMIN TAB TO TAB CONTAINER		
	this.tabContainer.addChild(this.settingsTab);
	this.tabContainer.selectChild(this.settingsTab);

	// SET ADD SOURCE ONCLICK
	dojo.connect(this.updateUserButton, "onClick", dojo.hitch(this, "updateUser"));

	// SET ADD SOURCE ONCLICK
	dojo.connect(this.addAwsButton, "onClick", dojo.hitch(this, "addAws"));

	//// ADD password CLASS TO PASSWORD INPUTS
	//dojo.addClass(this.newPassword.focusNode, "password");

	this.confirmPassword.parentWidget = this;
	this.newPassword.parentWidget = this;
	this.newPassword.target = this.confirmPassword;

	// SET DRAG SOURCE - LIST OF USERS
	this.initialiseSettings();
},

passwordsMatch : function () {
	var newPassword = this.newPassword.textbox.value;
	var confirmPassword = this.confirmPassword.textbox.value;

	return newPassword == confirmPassword;
},

addAws : function (event) {

	if ( this.savingAws == true )
	{
		return;
	}
	this.savingAws = true;

	// ALERT IF PASSWORDS DO NOT MATCH
	var amazonuserid = this.amazonuserid.value;
	amazonuserid = this.cleanEdges(amazonuserid);

	var volumesize = this.volumesize.value;
	if ( volumesize == null )
	{
		this.savingAws = false;
		return;
	}
	if ( ! volumesize.toString().match(/^\d+$/))
	{
		this.savingAws = false;
		return;
	}
	if ( volumesize < 10 || volumesize > 1000 )
	{
		this.savingAws = false;
		return;
	}

	if ( amazonuserid == "" )
	{
		this.savingAws = false;
		dojo.addClass(this.amazonuserid, 'invalid');
		return;
	}
	else
	{
		dojo.removeClass(this.amazonuserid, 'invalid');
	}


	// CLEAN UP WHITESPACE AND SUBSTITUTE NON-JSON SAFE CHARACTERS
	var query = new Object;
	query.username = Agua.cookie('username');
	query.volumesize = this.volumesize.value;
	query.amazonuserid = this.cleanEdges(this.amazonuserid.value);
	query.ec2publiccert = this.cleanEdges(this.ec2publiccert.value);
	query.ec2privatekey = this.cleanEdges(this.ec2privatekey.value);
	query.awsaccesskeyid = this.cleanEdges(this.awsaccesskeyid.value);
	query.secretaccesskey = this.cleanEdges(this.secretaccesskey.value);
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "addAws";

	var url = Agua.cgiUrl + "/agua?";

	// SEND TO SERVER
	dojo.xhrPut(
		{
			url: url,
			contentType: "json",
			putData: dojo.toJson(query),
			timeout: 15000,
			load: function(response, ioArgs) {
				return response;
			},
			error: function(response, ioArgs) {
				return response;
			}
		}
	);

	this.savingAws = false;
},

initialiseSettings : function () {

	// DISPLAY USERNAME
	var username = Agua.cookie('username');
	this.username.innerHTML = username;

	// INITIALISE USER INFO
	var user = Agua.getUser(username);
	this.firstname.value = user[1];
	this.lastname.value = user[2];
	this.email.value = user[3];

	// INITIALISE AWS SETTINGS
	var aws = Agua.getAws(username);

	this.amazonuserid.value = aws.amazonuserid || "";
	this.ec2publiccert.value = aws.ec2publiccert || "";
	this.ec2privatekey.value = aws.ec2privatekey || "";
	this.awsaccesskeyid.value = aws.awsaccesskeyid || "";
	this.secretaccesskey.value = aws.secretaccesskey || "";
},




// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS
reloadSettingsTabs : function () {

	var tabPaneNames = ["plugins.admin.GroupSettings"];
	for ( var i in tabPaneNames )
	{
		if ( this.parentWidget.paneWidgets[tabPaneNames[i]] != null )
		{
			this.parentWidget.paneWidgets[tabPaneNames[i]].reload();
		}
	}
},



setTable : function () {

	// DELETE EXISTING TABLE CONTENT
	while ( this.settingsTable.firstChild )
	{
		this.settingsTable.removeChild(this.settingsTable.firstChild);
	}
},




updateUser : function (event) {

	if ( this.savingUser == true )	return;
	this.savingUser = true;

	// GET PASSWORDS
	var oldPassword = this.cleanEdges(this.oldPassword.value);
	var newPassword = this.cleanEdges(this.newPassword.value);
	var confirmPassword = this.cleanEdges(this.confirmPassword.value);

	// RETURN IF ONLY SOME OF THE PASSWORDS ARE FILLED IN
	if ( ! (oldPassword  == '' && newPassword == '' && confirmPassword == '') )
	{
		if ( oldPassword  == '' ) { dojo.addClass(this.oldPassword, 'invalid'); }
		else { dojo.removeClass(this.oldPassword, 'invalid'); }
		if ( newPassword  == '' ) { dojo.addClass(this.newPassword, 'invalid'); }
		else { dojo.removeClass(this.newPassword, 'invalid'); }
		if ( confirmPassword  == '' ) { dojo.addClass(this.confirmPassword, 'invalid'); }
		else { dojo.removeClass(this.confirmPassword, 'invalid'); }		
		if ( oldPassword  == '' || newPassword == '' || confirmPassword == '' ) return;
	}

	// ALERT IF PASSWORDS DO NOT MATCH
	if ( newPassword != confirmPassword )
	{
		dojo.addClass(this.newPassword, 'invalid');
		dojo.addClass(this.confirmPassword, 'invalid');

		this.savingUser = false;
		return;
	}
	else
	{
		dojo.removeClass(this.newPassword, 'invalid');
		dojo.removeClass(this.confirmPassword, 'invalid');
	}

	// CLEAN UP WHITESPACE
	var user = new Object;
	user.username = Agua.cookie('username');
	user.firstname = this.cleanEdges(this.firstname.value);
	user.lastname = this.cleanEdges(this.lastname.value);
	user.email = this.cleanEdges(this.email.value);
	user.oldpassword = oldPassword;
	user.newpassword = newPassword;

	var url = Agua.cgiUrl + "/agua?";

	// CREATE JSON QUERY
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "updateUser";
	query.data = user;

	// SEND TO SERVER
	var thisObj = this;
	dojo.xhrPut(
		{
			url: url,
			contentType: "text",
			handleAs: "json",
			putData: dojo.toJson(query),
			timeout: 15000,
			load: function(response, ioArgs) {

				if ( response.error != null )
				{
					dojo.addClass(thisObj.oldPassword, 'invalid');
				}
				//return response;

			},
			error: function(response, ioArgs) {
				//return response;
			}
		}
	);

	this.savingUser = false;

}, // Settings.updateUser

// REMOVE WHITESPACE FROM EDGES OF TEXT
cleanEdges : function (string ) {
	if ( string == null )	{ 	return null; }
	string = string.replace(/^\s+/, '');
	string = string.replace(/\s+$/, '');
	return string;
}


}); // plugins.admin.Settings


