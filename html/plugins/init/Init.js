dojo.provide( "plugins.init.Init");

// INITIALISE AGUA - MOUNT DATA VOLUMES AND STORE ADMIN KEYS

// INHERITS
dojo.require("plugins.core.Common");
dojo.require("dijit.layout.ContentPane");
dojo.require("dijit.form.Button");
dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.Textarea");
dojo.require("dojo.parser");
dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dojo.parser");	// scan page for widgets and instantiate them

// FORM VALIDATION
dojo.require("plugins.form.ValidationTextarea");

dojo.declare(
    "plugins.init.Init",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "init/templates/init.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

//addingUser STATE
addingUser : false,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ "plugins/init/css/init.css" ],

// PARENT WIDGET
parentWidget : null,

// DEFAULT DATA VOLUME
defaultDataVolume : null,

/////}}

constructor : function(args) {
	// LOAD CSS
	this.loadCSS();		

	this.defaultDataVolume = args.dataVolume;
},

postCreate : function() {

	this.startup();
},


startup : function () {

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 

	this.datavolume.set("value", this.defaultDataVolume);

	dijit.Tooltip.defaultPosition = ['above', 'below'];

	// ADD ADMIN TAB TO TAB CONTAINER		
	dojo.byId("attachPoint").appendChild(this.initTab.domNode);

	// SET ADD SOURCE ONCLICK
	dojo.connect(this.saveUserButton, "onClick", dojo.hitch(this, "saveAws"));

	// SET ADD SOURCE ONCLICK
	dojo.connect(this.saveAwsButton, "onClick", dojo.hitch(this, "saveAws"));
},


toggleUserVolume : function (event) {

	if ( this.uservolumecheckbox.checked == true )
	{
		this.uservolume.set("disabled", false);
		dojo.addClass(this.uservolume.textbox, "enabled");
		this.uservolume.textbox.value = "";
		this.uservolume.textbox.focus();
	}
	if ( this.uservolumecheckbox.checked == false )
	{
		this.uservolume.set("disabled", true);
		this.uservolume.textbox.value = "  New volume";
		// WILL USE this.uservolume.value TO POPULATE QUERY
		this.uservolume.value = "";
		dojo.removeClass(this.uservolume.textbox, "enabled");
	}
},


toggleDataVolume : function (event) {

	if ( this.datavolumecheckbox.checked == true )
	{
		this.datavolume.set("disabled", false);
		dojo.addClass(this.datavolume.textbox, "enabled");
		this.datavolume.textbox.value = "";
		this.datavolume.textbox.focus();
	}
	if ( this.datavolumecheckbox.checked == false )
	{
		this.datavolume.set("disabled", true);
		this.datavolume.textbox.value = this.defaultDataVolume;
		// WILL USE this.datavolume.value TO POPULATE QUERY
		this.datavolume.value = this.defaultDataVolume;
		dojo.removeClass(this.datavolume.textbox, "enabled");
	}
},


saveAws : function (event) {

	if ( this.savingAws == true )
	{
		return;
	}
	this.savingAws = true;

	var parameters = [ "amazonuserid", "datavolume", "uservolume", "datavolumesize", "uservolumesize", "awsaccesskeyid", "secretaccesskey", "ec2publiccert", "ec2privatekey" ];
	for ( var i = 0; i < parameters.length; i++ )
	{
		var parameter = parameters[i];
		if ( this[parameter].isValid() == false )
		{
			if ( parameter == "datavolume" && this.datavolumecheckbox.checked == false )	continue;
			if ( parameter == "uservolume" && this.uservolumecheckbox.checked == false )	continue;
			this.savingAws = false;
			return;
		}
	}


	if ( this.uservolumecheckbox.checked == true )
	{
		if ( this.uservolume.isValid() == false )
		{
			this.savingAws = false;
			return;
		}
	}

	// CLEAN UP WHITESPACE AND SUBSTITUTE NON-JSON SAFE CHARACTERS
	var aws = new Object;
	aws.username = "admin";
	aws.amazonuserid = this.cleanEdges(this.amazonuserid.value);
	aws.datavolume = this.datavolume.value;
	aws.uservolume = this.uservolume.value;
	aws.datavolumesize = this.datavolumesize.value;
	aws.uservolumesize = this.uservolumesize.value;
	aws.datavolumecheckbox = this.datavolumecheckbox.checked == true ? 1 : 0;
	aws.uservolumecheckbox = this.uservolumecheckbox.checked == true ? 1 : 0;
	aws.ec2publiccert = this.cleanEdges(this.ec2publiccert.value);
	aws.ec2privatekey = this.cleanEdges(this.ec2privatekey.value);
	aws.awsaccesskeyid = this.cleanEdges(this.awsaccesskeyid.value);
	aws.secretaccesskey = this.cleanEdges(this.secretaccesskey.value);

	// ADD USER VOLUME IF 'CUSTOM' IS CHECKED
	if ( this.uservolumecheckbox.checked == true )
	{
		aws.uservolume = this.uservolume.value;
	}

	var url = this.cgiUrl + "/init.cgi?";

	// CREATE JSON QUERY
	var query = new Object;
	query.username = "agua";
	query.mode = "init";
	query.data = aws;

	// SEND TO SERVER
	var thisObj = this;
	dojo.xhrPut(
		{
			url: url,
			contentType: "text",
			putData: dojo.toJson(query),
			timeout: 15000,
			load: function(response, ioArgs) {

				thisObj.progress.innerHTML = response;
				dojo.addClass(thisObj.progress, "active");
				return response;
			},
			error: function(response, ioArgs) {
				return response;
			}
		}
	);

	this.savingAws = false;
},



// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS
reloadInitTabs : function () {

	var tabPaneNames = ["plugins.init.GroupInit"];
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
	while ( this.initTable.firstChild )
	{
		this.initTable.removeChild(this.initTable.firstChild);
	}
},


// REMOVE WHITESPACE FROM EDGES OF TEXT
cleanEdges : function (string) {
	if ( string == null )	{ 	return null; }
	string = string.replace(/^\s+/, '');
	string = string.replace(/\s+$/, '');
	return string;
}




}); // end of Init

