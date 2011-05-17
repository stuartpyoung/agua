dojo.provide("plugins.admin.Cluster");

// ALLOW THE ADMIN USER TO ADD, REMOVE AND MODIFY USERS

// NEW USERS MUST HAVE username AND email

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dijit.form.Button");
dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.Textarea");
dojo.require("dojo.parser");
dojo.require("dojo.dnd.Source");
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.admin.UserRow");

dojo.declare(
    "plugins.admin.Cluster",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/Cluster.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	//addingUser STATE
	addingUser : false,

	// OR USE @import IN HTML TEMPLATE
	cssFiles : [ "plugins/admin/css/Cluster.css" ],

	// PARENT WIDGET
	parentWidget : null,


	constructor : function(args)
	{
		// GET INFO FROM ARGS
		this.parentWidget = args.parentWidget;
		this.Cluster = args.parentWidget.Cluster;

        // LOAD SORIA AND FILEPICKER CSS
        this.loadCSS();		
	},

	postCreate : function()
	{

		this.startup();
	},


	startup : function ()
	{

		// COMPLETE CONSTRUCTION OF OBJECT
		this.inherited(arguments);	 

		// ADD ADMIN TAB TO TAB CONTAINER		
		this.tabContainer.addChild(this.ClusterTab);
		this.tabContainer.selectChild(this.ClusterTab);

		// SET ADD SOURCE ONCLICK
		dojo.connect(this.saveUserButton, "onClick", dojo.hitch(this, "saveUser"));

		// SET ADD SOURCE ONCLICK
		dojo.connect(this.saveAwsButton, "onClick", dojo.hitch(this, "saveAws"));

		// SET DRAG SOURCE - LIST OF USERS
		this.initialiseCluster();
	},


	saveAws : function (event)
	{

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
		var aws = new Object;
		aws.username = Agua.cookie('username');
		aws.volumesize = this.volumesize.value;
		aws.amazonuserid = this.cleanEdges(this.amazonuserid.value);
		aws.ec2publiccert = this.cleanEdges(this.ec2publiccert.value);
		aws.ec2privatekey = this.cleanEdges(this.ec2privatekey.value);
		aws.awsaccesskeyid = this.cleanEdges(this.awsaccesskeyid.value);
		aws.secretaccesskey = this.cleanEdges(this.secretaccesskey.value);

		var url = Agua.cgiUrl + "/agua?";

		// CREATE JSON QUERY
		var query = new Object;
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.mode = "saveAws";
		query.data = aws;

		// SEND TO SERVER
		dojo.xhrPut(
			{
				url: url,
				contentType: "text",
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

	initialiseCluster : function ()
	{

		// DISPLAY USERNAME
		var username = Agua.cookie('username');
		this.username.innerHTML = username;

		// INITIALISE USER INFO
		var user = Agua.getUser(username);
		this.firstname.value = user[1];
		this.lastname.value = user[2];
		this.email.value = user[3];

		// INITIALISE AWS Cluster
		var aws = Agua.getAws(username);

		this.amazonuserid.value = aws.amazonuserid || "";
		this.ec2publiccert.value = aws.ec2publiccert || "";
		this.ec2privatekey.value = aws.ec2privatekey || "";
		this.awsaccesskeyid.value = aws.awsaccesskeyid || "";
		this.secretaccesskey.value = aws.secretaccesskey || "";
	},




	// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS
	reloadClusterTabs : function ()
	{

		var tabPaneNames = ["plugins.admin.GroupCluster"];
		for ( var i in tabPaneNames )
		{
			if ( this.parentWidget.paneWidgets[tabPaneNames[i]] != null )
			{
				this.parentWidget.paneWidgets[tabPaneNames[i]].reload();
			}
		}
	},



	setTable : function ()
	{

		// DELETE EXISTING TABLE CONTENT
		while ( this.ClusterTable.firstChild )
		{
			this.ClusterTable.removeChild(this.ClusterTable.firstChild);
		}
	},




	saveUser : function (event)
	{

		if ( this.savingUser == true )
		{
			return;
		}
		this.savingUser = true;

		// ALERT IF PASSWORDS DO NOT MATCH
		var password = this.password.value;
		var passwordConfirm = this.passwordConfirm.value;
		password = this.cleanEdges(password);
		passwordConfirm = this.cleanEdges(passwordConfirm);
		if ( password != passwordConfirm )
		{
			dojo.addClass(this.password, 'invalid');
			dojo.addClass(this.passwordConfirm, 'invalid');

			this.savingUser = false;
			return;
		}
		else
		{
			dojo.removeClass(this.password, 'invalid');
			dojo.removeClass(this.passwordConfirm, 'invalid');
		}

		// CLEAN UP WHITESPACE AND SUBSTITUTE NON-JSON SAFE CHARACTERS
		var user = new Object;
		user.username = Agua.cookie('username');
		user.firstname = this.cleanEdges(this.firstname.value);
		user.email = this.cleanEdges(this.email.value);
		user.password = this.cleanEdges(this.password.value);

		var url = Agua.cgiUrl + "/agua?";

		// CREATE JSON QUERY
		var query = new Object;
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.mode = "saveUser";
		query.data = user;

		// SEND TO SERVER
		dojo.xhrPut(
			{
				url: url,
				contentType: "text",
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

		this.savingUser = false;

	}, // Cluster.saveUser

	// REMOVE WHITESPACE FROM EDGES OF TEXT
	cleanEdges : function (string)
	{
		if ( string == null )	{ 	return null; }
		string = string.replace(/^\s+/, '');
		string = string.replace(/\s+$/, '');
		return string;
	}




}); // plugins.admin.Cluster

