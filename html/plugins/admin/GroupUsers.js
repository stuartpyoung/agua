dojo.provide("plugins.admin.GroupUsers");

dojo.require("dijit.dijit"); // optimize: load dijit layer


// GENERAL FORM MODULES
dojo.require("dijit.form.Button");
dojo.require("dijit.form.ComboBox");
dojo.require("plugins.core.ComboBox");
dojo.require("dojo.data.ItemFileWriteStore");

// DRAG N DROP
dojo.require("dojo.dnd.Source");

// TOOLTIP
//dojo.require("dijit.Tooltip");

// SLIDER
dojo.require("dijit.form.Slider");

// PARSE
dojo.require("dojo.parser");

// INHERITS
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.admin.GroupUserRow");


dojo.declare(
    "plugins.admin.GroupUsers",

	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/groupusers.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	//addingSource STATE
	addingSource : false,

	// OR USE @import IN HTML TEMPLATE
	cssFiles : [ "plugins/admin/css/groupusers.css" ],

	// PARENT WIDGET
	parentWidget : null,

	// MAX. NO. OF USERS TO DISPLAY AT ANY ONE TIME IN DRAG SOURCE
	maxDisplayedUsers : 20,

 	// CONSTRUCTOR	
	constructor : function(args) {

		// REGISTER module path FOR PLUGINS


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
		this.tabContainer.addChild(this.groupusersTab);
		this.tabContainer.selectChild(this.groupusersTab);

		// SET GROUP COMBO
		this.setGroupCombo();

		// SET DRAG SOURCE - LIST OF SOURCES
		this.setDragSource();

		// SET DRAG SOURCE SLIDER
		this.setDragSourceSlider();

		// SET DROP TARGET - SOURCES ALREADY IN THE GROUP
		this.setDropTarget();

		// SET TRASH DROP TARGET
		this.setTrash();


		// example subscribe to events
		dojo.subscribe("/dnd/start", function(source){
			//console.debug("/dnd/start Starting the drop", source);
		});
		dojo.subscribe("/dnd/drop/before", function(source, nodes, copy, target){
			//console.debug(copy ? "/dnd/drop/before Copying from" : "Moving from", source, "to", target, "before", target.before);
		});
		dojo.subscribe("/dnd/drop", function(source, nodes, copy, target){
			//console.debug(copy ? "/dnd/drop Copying from" : "Moving from", source, "to", target, "before", target.before);
		});

	},


	// reload
	//
	// RELOAD THE GROUP COMBO AND DRAG SOURCE
	// (CALLED AFTER CHANGES TO SOURCES OR GROUPS DATA IN OTHER TABS)
	//
	reload : function ()
	{
		// SET GROUP COMBO
		this.setGroupCombo();

		// SET DRAG SOURCE - LIST OF SOURCES
		this.setDragSource();

		// SET DROP TARGET - SOURCES ALREADY IN THE GROUP
		this.setDropTarget();
	},



	// SET GROUP COMBO BOX
	setGroupCombo : function ()
	{
		// GET GROUP NAMES		
		var groupNames = Agua.getGroupNames();

		// SET STORE
		var data = {identifier: "name", items: []};
		for ( var i = 0; i < groupNames.length; i++ )
		{
			data.items[i] = { name: groupNames[i]	};
		}
		var store = new dojo.data.ItemFileWriteStore(	{	data: data	}	);

		this.groupCombo.popupClass = "groupusers groupCombo dijitReset dijitMenu";
		this.groupCombo.wrapperClass = "groupusers dijitPopup";
		this.groupCombo.itemHeight = 30;

		// SET COMBO
		this.groupCombo.store = store;
		this.groupCombo.startup();

		// SET COMBO VALUE
		var firstValue = groupNames[0];
		this.groupCombo.setValue(firstValue);

		// CONNECT ONCLICK WITH dojo.connect TO BUILD TABLE
		var groupUsersObjects = this;
		dojo.connect(this.groupCombo, "onChange", function(event) {
			groupUsersObjects.setDropTarget();
		});
	},



	setDragSource : function (position)
	{

		// REMOVE ALL EXISTING CONTENT
		while ( this.dragSourceContainer.firstChild )
		{
			if ( dijit.byNode(this.dragSourceContainer)
				&& dijit.byNode(this.dragSourceContainer).destroy )
			{
				dijit.byNode(this.dragSourceContainer).destroy();
			}
			this.dragSourceContainer.removeChild(this.dragSourceContainer.firstChild);
		}

		// SET position IF NOT DEFINED
		if ( position == null )	position = 0;

		var dataArray = new Array;
		var userArray = Agua.getUsers();

		// RETURN IF USER ARRAY IS NULL OR EMPTY
		if ( userArray == null || userArray.length == 0 )
		{
			return;
		}

		var MAXUSERS = this.maxDisplayedUsers;
		var MULTIPLE = ( userArray.length - MAXUSERS ) / 100;
		MULTIPLE = MULTIPLE ? MULTIPLE : 1;
		var start = parseInt(position * MULTIPLE);
		var end = parseInt( (position * MULTIPLE) + MAXUSERS );
		if ( ! end || end > userArray.length )	end = userArray.length;

		// GENERATE USER DATA TO INSERT INTO DND USER TABLE
		for ( var j = start; j < end; j++ )
		{
			var data = userArray[j];				
			data.toString = function () { return this.name; }
			dataArray.push( { data: data, type: ["draggableUser"] } );
		}

		// GENERATE DND SOURCE
		var dragSource = new dojo.dnd.Source(
			this.dragSourceContainer,
			{
				copyOnly: true,
				selfAccept: false,
				accept : [ "none" ]
			}
		);
		dragSource.insertNodes(false, dataArray);

		// users TABLE QUERY:
		// SELECT DISTINCT username, firstname, lastname, email, description FROM users ORDER BY username};

		// SET TABLE ROW STYLE IN dojDndItems
		var allNodes = dragSource.getAllNodes();
		for ( var k = 0; k < allNodes.length; k++ )
		{
			// ADD CLASS FROM type TO NODE
			var node = allNodes[k];

			// SET NODE name AND description
			node.username = dataArray[k].data[0] ? dataArray[k].data[0]: '';
			node.firstname = dataArray[k].data[1] ? dataArray[k].data[1]: '';
			node.lastname = dataArray[k].data[2] ? dataArray[k].data[2]: '';
			node.email = dataArray[k].data[3] ? dataArray[k].data[3]: '';
			node.description = dataArray[k].data[4] ? dataArray[k].data[4]: '';

			node.firstname = this.firstLetterUpperCase(node.firstname);
			node.lastname = this.firstLetterUpperCase(node.lastname);

			// CHECK ATTRIBUTES ARE NOT NULL
			if ( node.description == null )
				node.description = '';

			if ( node.username != null && node.username != '' )
			{
				var user = {
					username : node.username,
					email : node.email,
					firstname : node.firstname,
					lastname : node.lastname,
					description : node.description,
					location : ''
					//datetime : node.datetime,
				};
				user.parentWidget = this;			

				var groupUserRow = new plugins.admin.GroupUserRow(user);
				node.innerHTML = '';
				node.appendChild(groupUserRow.domNode);
			}
		}

		var groupObject = this;
		dragSource.creator = function (item, hint)
		{

			var node = dojo.doc.createElement("div");
			node.username = item[0];
			node.firstname = item[1];
			node.lastname = item[2];
			node.email = item[3];
			node.description = item[4];

			node.id = dojo.dnd.getUniqueId();
			node.className = "dojoDndItem";


			// SET FANCY FORMAT IN NODE INNERHTML
			node.innerHTML = "<table> <tr><td><strong style='color: darkred'>" + item[0] + "</strong></td></tr><tr><td> " + item[3] + "</td></tr></table>";

			return {node: node, data: item, type: ["text"]};
		};
	},


	setDragSourceSlider : function ()
	{

		// ONMOUSEUP
		var groupUsersObject = this;
		dojo.connect(this.dragSourceSlider, "onMouseUp", dojo.hitch(this, function(e)
		{
			var position = parseInt(this.dragSourceSlider.getValue());
			this.setDragSource(position);
		}));
	},


	setDropTarget : function ()
	{

		// DELETE EXISTING CONTENTS OF DROP TARGET
		while ( this.dropTargetContainer.firstChild )
		{
			this.dropTargetContainer.removeChild(this.dropTargetContainer.firstChild);
		}


		//// GET THE SOURCES IN THIS GROUP
		var groupname = this.groupCombo.getValue();
		var sourceArray = Agua.getGroupUsers();
		sourceArray = this.filterByKeyValues(sourceArray, ["groupname"], [groupname]);

		if ( sourceArray == null )
		{
			return;
		}		

		var dataArray = new Array;
		for ( var j = 0; j < sourceArray.length; j++ )
		{
			var data = sourceArray[j];				
			data.toString = function () { return this.name; }
			var description = sourceArray[j].description;

			dataArray.push( { data: data, type: ["draggableUser"], description: description  } );
		}

		// GENERATE DROP TARGET
		var dropTarget = new dojo.dnd.Source(
			this.dropTargetContainer,
			{
				//copyOnly: true,
				//selfAccept: false
				//,
				accept : [ "draggableUser" ]
			}
		);
		dropTarget.insertNodes(false, dataArray );

		// SET TABLE ROW STYLE IN dojDndItems
		var allNodes = dropTarget.getAllNodes();
		for ( var k = 0; k < allNodes.length; k++ )
		{


			// ADD CLASS FROM type TO NODE
			var node = allNodes[k];
			var nodeClass = dataArray[k].type;
			dojo.addClass(node, nodeClass);

			// SET NODE ATTRIBUTES
			node.name = dataArray[k].data.name;
			node.description = dataArray[k].data.description;
			node.location = dataArray[k].data.location;

			var user = {
				username : node.name,
				email : '',
				location : node.location,
				firstname : '',
				lastname : '',
				description : node.description
			};
			user.parentWidget = this;			

			var groupUserRow = new plugins.admin.GroupUserRow(user);
			node.innerHTML = '';
			node.appendChild(groupUserRow.domNode);
		}

		var groupUsersObject = this;
		dropTarget.creator = function (item, hint)
		{

			var username;
			var description;
			var email;


			// IF item IS AN ARRAY, IT CAME FROM THE DRAG SOURCE
			if ( item.length != null )
			{
				username = item[0];
				description = groupUsersObject.firstLetterUpperCase(item[1]) + "  " + groupUsersObject.firstLetterUpperCase(item[2]);
				email = item[3];
			}
			// OTHERWISE, IT CAME FROM THE DROP TARGET ITSELF
			else
			{
				username = item.name;
				description = item.description;
				email = item.location;
			}

			// ADD VALUES TO NODE SO THAT THEY GET PASSED TO
			// this.addToGroup
			var node = dojo.doc.createElement("div");
			node.name = username;
			node.description = description;
			node.email = email;
			node.location = '';
			node.id = dojo.dnd.getUniqueId();
			node.className = "dojoDndItem";


			var user = {
				username : node.name,
				email : node.email,
				location : node.location,
				firstname : '',
				lastname : '',
				description : node.description
			};

			user.parentWidget = this;			

			var groupUserRow = new plugins.admin.GroupUserRow(user);
			node.innerHTML = '';
			node.appendChild(groupUserRow.domNode);


			// ACCEPTABLE NODE DROPPED ONTO SELF --> ADD TO USER GROUP.
			if ( hint != 'avatar' )
			{
				groupUsersObject.addToGroup(node.name, node.description, node.email);
			}

			return {node: node, data: item, type: ["draggableUser"]};
		};


		// ADD NODE IF DROPPED FROM OTHER SOURCE, DELETE IF DROPPED FROM SELF
		dojo.connect(dropTarget, "onDndDrop", function(source, nodes, copy, target){


			// NODE DROPPED FROM SELF --> DELETE NODE
			if ( source == this && target != this )
			{
				for ( var i = 0; i < nodes.length; i++ )
				{
					var node = nodes[i];


					groupUsersObject.removeFromGroup(node.name, node.description, node.location);



				}
			}

			// DO NOTHING IF NODE WAS DROPPED FROM DRAG SOURCE
			// ( IT HAS ALREADY BEEN GENERATED BY dropTarget.creator() )
			else
			{
				for ( var i = 0; i < nodes.length; i++ )
				{
					var node = nodes[i];
				}

			}


			// REMOVE DUPLICATE NODES
			var currentNodes = dropTarget.getAllNodes();

			////console.dir(currentNodes);
			if ( currentNodes == null || ! currentNodes )
			{
				return;
			}

			// DELETE DUPLICATE NODES
			var names = new Object;
			for ( var i = 0; i < currentNodes.length; i++ )
			{
				var node = currentNodes[i];
				if ( ! names[node.name] )
				{
					names[node.name] = 1;
				}
				else
				{
					//document.removeChild(node); // Node was not found" code: "8

					// HACK TO AVOID THIS ERROR: node.parentNode is null
					try {
						//i--;
						node.parentNode.removeChild(node)
					}
					catch (e) {
					}

//						currentNodes.removeChild(node); //currentNodes.removeChild is not a function
					//currentNodes.splice(i, 1); // no effect on nodes in dnd.Source
				}

			}
		});
	},

	addToGroup : function (name, description, location )
	{

		var groupObject = new Object;
		groupObject.username = Agua.cookie('username');
		groupObject.name = name;
		groupObject.description = description;
		groupObject.location = location;

		// ADD SOURCE OBJECT TO THE SOURCES IN THIS GROUP
		var groupName = this.groupCombo.getValue();

		if ( Agua.addUserToGroup(groupName, groupObject) == false )
		{
			return;
		}

		// ADD THE SOURCE INTO THE groupusers TABLE ON THE SERVER
		var data = new Object;
		data.name = name;
		data.description = description;
		data.location = location;
		data.groupname = groupName;
		data.type = "user";

		var url = Agua.cgiUrl + "agua?";
		var query = new Object;
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.data = data;
		query.mode = "addToGroup";

		var queryString = dojo.toJson(query).replace(/undefined/g, '""');

		// SEND TO SERVER
		dojo.xhrPut(
			{
				url: url,
				contentType: "text",
				sync : false,
				handleAs: "json",
				putData: queryString,
				timeout: 15000,
				load: function(data)
				{
				},
				error: function(response, ioArgs) {
					return response;
				}
			}
		);
	},



	removeFromGroup : function (name, description, location )
	{

		// REMOVE SOURCE FROM THE SOURCES IN THIS GROUP
		var groupName = this.groupCombo.getValue();

		var groupObject = new Object;
		groupObject.name = name;
		groupObject.description = description;
		groupObject.location = location;

		if ( Agua.removeUserFromGroup(groupName, groupObject) == false )
		{
			return;
		}

		// REMOVE THE SOURCE FROM THE groupusers TABLE ON THE SERVER
		var data = new Object;
		data.name = String(name);
		data.description = description;
		data.location = location;
		data.groupname = groupName;
		data.type = "user";

		var url = Agua.cgiUrl + "agua?";
		var query = new Object;
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.data = data;
		query.mode = "removeFromGroup";

		var queryString = dojo.toJson(query).replace(/undefined/g, '""');

		// SEND TO SERVER
		dojo.xhrPut(
			{
				url: url,
				contentType: "text",
				sync : false,
				handleAs: "json",
				putData: queryString,
				timeout: 15000,
				load: function(data)
				{
				},
				error: function(response, ioArgs) {
					return response;
				}
			}
		);
	},


	//	DELETE NODE IF DROPPED INTO TRASH. 
	//	REMOVAL OF DATA IS ACCOMPLISHED IN THE
	// onDndDrop LISTENER OF THE SOURCE
	setTrash : function ()
	{

		var trash = new dojo.dnd.Source(
			this.trashContainer,
			{
				accept : [ "draggableUser" ]
			}
		);

		// REMOVE DUPLICATE NODES
		dojo.connect(trash, "onDndDrop", function(source, nodes, copy, target){


			// NODE DROPPED ON SELF --> DELETE THE NODE
			if ( target == this )
			{
				var currentNodes = trash.getAllNodes();
				for ( var i = 0; i < currentNodes.length; i++ )
				{
					var node = currentNodes[i];
					// HACK TO AVOID THIS ERROR: node.parentNode is null
					try {
						node.parentNode.removeChild(node)
					}
					catch (e) {
					}
				}
			}
		});
	}

}); // plugins.admin.GroupUsers
