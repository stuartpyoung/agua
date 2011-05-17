dojo.provide("plugins.admin.Users");

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
    "plugins.admin.Users",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "admin/templates/users.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

//addingUser STATE
addingUser : false,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ "plugins/admin/css/users.css" ],

// PARENT WIDGET
parentWidget : null,

/////}}}

constructor : function(args) {
	// GET INFO FROM ARGS
	this.parentWidget = args.parentWidget;
	this.users = args.parentWidget.users;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();		
},

postCreate : function() {

	this.startup();
},


startup : function() {

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 

	// ADD ADMIN TAB TO TAB CONTAINER		
	this.tabContainer.addChild(this.usersTab);
	this.tabContainer.selectChild(this.usersTab);

	// SET DRAG SOURCE - LIST OF USERS
	this.setTable();

	// SET NEW SOURCE FORM
	this.setNewUser();

	// SET TRASH DROP TARGET
	this.setTrash();	
},



// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS
reloadUserTabs : function () {

	var tabPaneNames = ["plugins.admin.GroupUsers"];
	for ( var i in tabPaneNames )
	{
		if ( this.parentWidget.paneWidgets[tabPaneNames[i]] != null )
		{
			this.parentWidget.paneWidgets[tabPaneNames[i]].reload();
		}
	}
},


clearValue : function (textarea) {

	if ( textarea.clicked == true ) return;

	textarea.clicked = true;
	textarea.value = '';
	textarea.focus();
},

setNewUser : function () {

	// FOCUS ON WIDGET TO FIX ITS WIDTH
	this.newUsername.focus();

	// SET ADD SOURCE ONCLICK
	dojo.connect(this.addUserButton, "onClick", dojo.hitch(this, "addUser"));

	// SET ONCLICK TO CANCEL DEFAULT TEXT
	dojo.connect(this.newUsername, "onclick", dojo.hitch(this, "clearValue", this.newUsername));
	dojo.connect(this.newUsername, "onfocus", dojo.hitch(this, "clearValue", this.newUsername));
	dojo.connect(this.newFirstname, "onclick", dojo.hitch(this, "clearValue", this.newFirstname));
	dojo.connect(this.newFirstname, "onfocus", dojo.hitch(this, "clearValue", this.newFirstname));
	dojo.connect(this.newLastname, "onclick", dojo.hitch(this, "clearValue", this.newLastname));
	dojo.connect(this.newLastname, "onfocus", dojo.hitch(this, "clearValue", this.newLastname));
	dojo.connect(this.newEmail, "onclick", dojo.hitch(this, "clearValue", this.newEmail));
	dojo.connect(this.newEmail, "onfocus", dojo.hitch(this, "clearValue", this.newEmail));		
	dojo.connect(this.newPassword, "onclick", dojo.hitch(this, "clearValue", this.newPassword));
	dojo.connect(this.newPassword, "onfocus", dojo.hitch(this, "clearValue", this.newPassword));		
},

setTable : function () {

	// DELETE EXISTING TABLE CONTENT
	while ( this.usersTable.firstChild )
	{
		this.usersTable.removeChild(this.usersTable.firstChild);
	}

	var dataArray = new Array;
	var userArray = Agua.getUsers();

	userArray = this.sortHasharray(userArray, 'username');

	// CHECK userArray IS NOT NULL OR EMPTY
	if ( userArray == null || userArray.length == 0 )
	{
		return;
	}

	// GENERATE dataArray TO INSERT INTO DND SOURCE TABLE
	for ( var j = 0; j < userArray.length; j++ )
	{
		var data = new Object;
		data.username = userArray[j][0];				
		data.firstname = userArray[j][1];				
		data.lastname = userArray[j][2];				
		data.email = userArray[j][3];				
		data.password = "its-a-secret";
		data.toString = function () { return this.username; }
		dataArray.push( { data: data, type: ["draggableItem"] } );
	}

	// GENERATE DND SOURCE
	var dragSource = new dojo.dnd.Source(
		this.usersTable,
		{
			copyOnly: true,
			selfAccept: false,
			accept : [ "none" ]
		}
	);
	dragSource.insertNodes(false, dataArray);

	// SET TABLE ROW STYLE IN dojDndItems
	var allNodes = dragSource.getAllNodes();

	for ( var k = 0; k < allNodes.length; k++ )
	{
		// ADD CLASS FROM type TO NODE
		var node = allNodes[k];


		// SET NODE username AND firstname
		node.username = dataArray[k].data.username;
		node.firstname = dataArray[k].data.firstname;
		node.lastname = dataArray[k].data.lastname;
		node.email = dataArray[k].data.email;
		node.password = dataArray[k].data.password;
		if ( node.firstname == null )	{	node.firstname = '';	}
		if ( node.lastname == null )	{	node.lastname = '';	}
		if ( node.email == null )	{	node.email = '';	}

		var userRow = new plugins.admin.UserRow({
			username : node.username,
			firstname : node.firstname,
			lastname : node.lastname,
			email : node.email,
			password : node.password,
			parentWidget : this
		});

		node.innerHTML = '';
		node.appendChild(userRow.domNode);

	}

	var userObject = this;
	//dojo.connect(dragSource, "creator", userObject.specialAvatar );

	dragSource.creator = function (item, hint)
	{

		var node = dojo.doc.createElement("div");
		node.username = item.username;
		node.firstname = item.firstname;
		node.lastname = item.lastname;
		node.email = item.email;
		node.password = item.password;
		node.id = dojo.dnd.getUniqueId();
		node.className = "dojoDndItem";


		// SET FANCY FORMAT IN NODE INNERHTML
		node.innerHTML = "<table> <tr><td colspan=2><strong style='color: darkred'>" + item.username + "</strong></td></tr><tr><td> " + item.firstname + "</td><td> " + item.lastname + "</td></tr></table>";

		return {node: node, data: item, type: ["text"]};
	};


},


editUserRow : function (userRow, node) {

	// RETURN IF ALREADY EDITING SOURCE ROW (I.E., MULTIPLE CLICKS)
	if ( this.editingUserRow == true ) return;
	this.editingUserRow = true;

	// REPLACE THE TD INNERHTML WITH A TEXTAREA
	var text = node.innerHTML;

	node.innerHTML = '';
	if ( text == null || ! text ) text = '';

	// RETURN IF THIS IS A DOUBLE-CLICK
	if ( text.match(/^<i/) ||
		text.match(/^<br/) ||
		text.match(/^<fieldset/) ||
		text.match(/^<textarea/) )
	{
		this.editingUserRow = false;
		return;
	}

	// CREATE INPUT TEXT AREA
	var textarea = document.createElement('textarea');
	dojo.addClass(textarea, 'editUserRow');
	node.appendChild(textarea);
	textarea.value = text;
	textarea.focus();

	// SET NEW PROJECT LISTENER
	var thisObject = this;
	dojo.connect(textarea, "onkeypress", function(evt){

		// summary: handles keyboard events
		var key = evt.charOrCode;

		if ( key == 13 )
		{
			var newText = textarea.value;

			var user = new Object;
			user.username	=	userRow.username.innerHTML;
			user.firstname	=	userRow.firstname.innerHTML;
			user.lastname	=	userRow.lastname.innerHTML;
			user.email		=	userRow.email.innerHTML;
			user.password	=	userRow.password.innerHTML;

			// AVOID NESTED CLICKS
			if ( user.firstname.match(/^<textarea/) )
				user.firstname = userRow.firstname.firstChild.value;

			if ( user.lastname.match(/^<textarea/) )
				user.lastname = userRow.lastname.lastChild.value;

			if ( user.email.match(/^<textarea/) )
				user.email = userRow.email.firstChild.value;

			if ( user.password.match(/^<textarea/) )
				user.password = userRow.password.firstChild.value;

			// REMOVE WHITESPACE
			//user.username = user.username.replace(/^\s+/, '');
			//user.username = user.username.replace(/\s+$/, '');
			user.username = user.username.match(/^\s*(.+?)\s*$/)[1];

			user.firstname = user.firstname.replace(/^\s+/, '');
			user.firstname = user.firstname.replace(/\s+$/, '');
			user.lastname = user.lastname.replace(/^\s+/, '');
			user.lastname = user.lastname.replace(/\s+$/, '');
			user.email = user.email.replace(/^\s+/, '');
			user.email = user.email.replace(/\s+$/, '');
			user.password = user.password.replace(/^\s+/, '');
			user.password = user.password.replace(/\s+$/, '');

			if ( user.username != '' )
				//&& user.firstname != '' 
				//&& user.email != '' )
			{
				// REMOVE ORIGINAL user OBJECT FROM Agua.users ARRAY
				// THEN ADD NEW SOURCE OBJECT TO Agua.users ARRAY
				Agua.removeUser({ username: user.username});
				Agua.addUser(user);

				// REMOVE TEXTAREA
				node.removeChild(textarea);
				node.innerHTML = newText;

				// CANCEL EDITING
				thisObject.editingUserRow = false;

				// REDO SOURCE TABLE
				//thisObject.setTable();

				// SAVE NEW SOURCE TO REMOTE DATABASE
				thisObject.saveUser(user);
			}
		}
		else if (key == dojo.keys.ESCAPE)
		{
			projectsObject.editingProjectRow = false;

			// REMOVE TEXTAREA
			node.removeChild(textarea);
			if ( text == null || text == '' )
			{
				text = "<i class='default'>Project notes</i>";
			}
			node.innerHTML = text;

		}
	});

},


// REMOVE USER FROM this.users
removeUser : function (username) {

	// CLEAN UP WHITESPACE
	username = username.replace(/\s+$/,'');
	username = username.replace(/^\s+/,'');

	var userObject = { username: username };

	// REMOVING SOURCE FROM Agua.users
	var success = Agua.removeUser(userObject)

	// RESET THE USERS TABLE
	this.setTable();

	var url = Agua.cgiUrl + "/agua?";

	// CREATE JSON QUERY
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "removeUser";
	query.data = userObject;

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

	// RELOAD SOURCE-RELATED TABS
	this.reloadUserTabs();

}, // Users.removeUser


addUser : function (node) {

	var username = this.newUsername.value;
	var firstname = this.newFirstname.value;
	var lastname = this.newLastname.value;
	var email = this.newEmail.value;
	var password = this.newPassword.value;

	// CHECK FOR VALID INPUTS
	if ( username == '' || username.match(/^\s*Username\s*$/) )
	{
		dojo.addClass(this.newUsername, 'invalid');
	}
	else{
		dojo.removeClass(this.newUsername, 'invalid');
	}
	if ( firstname == '' || firstname.match(/^\s*Firstname\s*$/) )
	{
		dojo.addClass(this.newFirstname, 'invalid');
	}
	else{
		dojo.removeClass(this.newFirstname, 'invalid');
	}
	if ( firstname == '' || firstname.match(/^\s*Lastname\s*$/) )
	{
		dojo.addClass(this.newLastname, 'invalid');
	}
	else{
		dojo.removeClass(this.newLastname, 'invalid');
	}

	if ( email == '' || email.match(/^\s*Email\s*$/) )
	{
		dojo.addClass(this.newEmail, 'invalid');
	}
	else{
		dojo.removeClass(this.newEmail, 'invalid');
	}

	if ( password == '' || password.match(/^\s*Password\s*$/) )
	{
		dojo.addClass(this.newPassword, 'invalid');
	}
	else{
		dojo.removeClass(this.newPassword, 'invalid');
	}

	if ( username == '' || username.match(/^\s*Name\s*$/) 
		|| firstname == '' || firstname.match(/^\s*Firstname\s*$/)
		|| lastname == '' || lastname.match(/^\s*Lastname\s*$/)
		|| email == '' || email.match(/^\s*Email\s*$/) 
		|| password == '' || password.match(/^\s*Password\s*$/) )
	{
		return;
	}

	var userObject = {
		username: username,
		firstname: firstname,
		lastname: lastname,
		email: email,
		password: password
	};

	if ( Agua.isUser(userObject) )
	{
		return;
	}

	// ADD USER TO USERS ARRAY
	var added = Agua.addUser(userObject);
	if ( added == false )
	{
		return false;
	}

	// RESET SOURCE TABLE
	this.setTable();

	if ( this.savingUser == true )
	{
		return;
	}
	this.savingUser = true;

	// CLEAN UP WHITESPACE AND SUBSTITUTE NON-JSON SAFE CHARACTERS
	userObject.originalName = this.jsonSafe(node.originalName, 'toJson');
	userObject.username = this.jsonSafe(userObject.username, 'toJson');
	userObject.firstname = this.jsonSafe(userObject.firstname, 'toJson');
	userObject.email = this.jsonSafe(userObject.email, 'toJson');
	userObject.password = this.jsonSafe(userObject.password, 'toJson');

	var url = Agua.cgiUrl + "/agua?";

	// CREATE JSON QUERY
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "addUser";
	query.data = userObject;
	console.dir(query);

	// SEND TO SERVER
	var thisObj = this;
	dojo.xhrPut(
		{
			url: url,
			contentType: "text",
			putData: dojo.toJson(query),
			timeout: 15000,
			load: function(response, ioArgs) {
				// RELOAD SOURCE-RELATED TABS
				thisObj.reloadUserTabs();

				return response;
			},
			error: function(response, ioArgs) {

				// RELOAD SOURCE-RELATED TABS
				thisObj.reloadUserTabs();

				return response;
			}
		}
	);

	this.savingUser = false;

},	// Users.addUser



saveUser : function (user)  {

	if ( this.savingUser == true )
	{
		return;
	}
	this.savingUser = true;

	// CLEAN UP WHITESPACE AND SUBSTITUTE NON-JSON SAFE CHARACTERS
	user.originalName = this.jsonSafe(user.originalName, 'toJson');
	user.username = this.jsonSafe(user.username, 'toJson');
	user.firstname = this.jsonSafe(user.firstname, 'toJson');
	user.email = this.jsonSafe(user.email, 'toJson');
	user.password = this.jsonSafe(user.password, 'toJson');

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

}, // Users.saveUser




// setTrash
//
//	DELETE NODE IF DROPPED INTO TRASH. ACTUAL REMOVAL FROM THE
//	DATA IS ACCOMPLISHED IN THE onDndDrop LISTENER OF THE SOURCE
//
setTrash : function () {

	var trash = new dojo.dnd.Source(
		this.trashContainer,
		{
			accept : [ "draggableItem" ]
		}
	);

	// REMOVE DUPLICATE NODES
	var thisObject = this;
	dojo.connect(trash, "onDndDrop", function(user, nodes, copy, target){
		// NODE DROPPED ON SELF --> DELETE THE NODE
		if ( target == this )
		{

			for ( var i = 0; i < nodes.length; i++ )
			{
				var node = nodes[i];

				// HACK TO AVOID THIS ERROR: node.parentNode is null
				try {
					node.parentNode.removeChild(node);

					thisObject.removeUser(node.username);
				}
				catch (e) {
				}
			}

			// DELETE EXISTING TABLE CONTENT
			//console.dir(thisObject.trashContainer);
			while ( thisObject.trashContainer.childNodes.length > 2 )
			{

				thisObject.trashContainer.removeChild(thisObject.trashContainer.childNodes[2]);
			}


		}
	});
}



}); // plugins.admin.Users

