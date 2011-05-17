dojo.provide("plugins.admin.GroupProjects");


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
dojo.require("plugins.admin.GroupProjectRow");


dojo.declare(
    "plugins.admin.GroupProjects",

	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/groupprojects.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	//addingSource STATE
	addingSource : false,

	// OR USE @import IN HTML TEMPLATE
	cssFiles : [ "plugins/admin/css/groupprojects.css" ],

	// PARENT WIDGET
	parentWidget : null,

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
		this.tabContainer.addChild(this.groupprojectsTab);
		this.tabContainer.selectChild(this.groupprojectsTab);

		// SET GROUP COMBO
		this.setGroupCombo();

		// SET SOURCE MANAGER BUTTON
		//this.sourceManagerButton();

		// SET DRAG SOURCE - LIST OF SOURCES
		this.setDragSource();

		// SET DRAG SOURCE SLIDER
		this.setDragSourceSlider();

		// SET DROP TARGET - SOURCES ALREADY IN THE GROUP
		this.setDropTarget();

		// SET TRASH DROP TARGET
		this.setTrash();


		//// example subscribe to events
		//dojo.subscribe("/dnd/start", function(source){
		//	//console.debug("/dnd/start Starting the drop", source);
		//});
		//dojo.subscribe("/dnd/drop/before", function(source, nodes, copy, target){
		//	//console.debug(copy ? "/dnd/drop/before Copying from" : "Moving from", source, "to", target, "before", target.before);
		//});
		//dojo.subscribe("/dnd/drop", function(source, nodes, copy, target){
		//	//console.debug(copy ? "/dnd/drop Copying from" : "Moving from", source, "to", target, "before", target.before);
		//});
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

		this.groupCombo.popupClass = "groupprojects groupCombo dijitReset dijitMenu";
		this.groupCombo.wrapperClass = "groupprojects dijitPopup";
		this.groupCombo.itemHeight = 30;

		// SET COMBO
		this.groupCombo.store = store;
		this.groupCombo.startup();

		// SET COMBO VALUE
		var firstValue = groupNames[0];
		this.groupCombo.setValue(firstValue);

		// CONNECT ONCLICK WITH dojo.connect TO BUILD TABLE
		var groupProjectsObjects = this;
		dojo.connect(this.groupCombo, "onChange", function(event) {
			groupProjectsObjects.setDropTarget();
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

		// RETURN IF USER ARRAY IS NULL OR EMPTY
		if ( Agua.projects == null || Agua.projects.length == 0 )
		{
			return;
		}

		var projectArray = Agua.projects;

		// GENERATE USER DATA TO INSERT INTO DND USER TABLE
		var MAXUSERS = this.maxDisplayedProjects;
		var MULTIPLE = ( projectArray.length - MAXUSERS ) / 100;
				MULTIPLE = MULTIPLE ? MULTIPLE : 1;

		var start = parseInt(position * MULTIPLE);
		var end = parseInt( (position * MULTIPLE) + MAXUSERS );
		if ( ! end || end > projectArray.length )	end = projectArray.length;


		// SORT PROJECT ARRAY

		// GENERATE USER DATA TO INSERT INTO DND USER TABLE
		var dataArray = new Array;
		for ( var j = start; j < end; j++ )
		{
			var data = projectArray[j];				
			data.toString = function () { return this.name; }
			dataArray.push( { data: data, type: ["draggableItem"] } );
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

		// projects TABLE QUERY:
		// SELECT DISTINCT projectname, firstname, lastname, email, description FROM projects ORDER BY projectname};
		// SET TABLE ROW STYLE IN dojDndItems
		var allNodes = dragSource.getAllNodes();
		for ( var k = 0; k < allNodes.length; k++ )
		{

			// ADD CLASS FROM type TO NODE
			var node = allNodes[k];

			// SET NODE name AND description
			node.name = dataArray[k].data.name;
			node.description = dataArray[k].data.description;
			if ( node.description == null ) node.description = '';

			// CHECK ATTRIBUTES ARE NOT NULL
			if ( node.description == null )
				node.description = '';

			if ( node.name != null && node.name != '' )
			{
				var project = {
					name : node.name,
					description : node.description
				};
				project.parentWidget = this;			

				var groupProjectRow = new plugins.admin.GroupProjectRow(project);
				//groupProjectRow.toggle();
				node.innerHTML = '';
				node.appendChild(groupProjectRow.domNode);
			}
		}

		var sourceObject = this;
		dragSource.creator = function (item, hint)
		{

			var node = dojo.doc.createElement("div");
			node.name = item.name;
			node.description = item.description;

			node.id = dojo.dnd.getUniqueId();
			node.className = "dojoDndItem";


			// SET FANCY FORMAT IN NODE INNERHTML
			node.innerHTML = "<table> <tr><td> <img src='http://localhost/agua/0.4/plugins/admin/images/users-18.png' <strong style='color: darkred'>" + node.name + "</strong></td></tr><tr><td> " + node.description + "</td></tr></table>";

			return {node: node, data: item, type: ["text"]};
		};
	},


	setDragSourceSlider : function ()
	{

		// ONMOUSEUP
		var groupProjectsObject = this;
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
		var projectArray = Agua.getProjectsByGroup(groupname);
		if ( projectArray == null )
		{
			return;
		}		

		var dataArray = new Array;
		for ( var j = 0; j < projectArray.length; j++ )
		{
			var data = projectArray[j];				
			data.toString = function () { return this.name; }
			var description = projectArray[j].description;

			dataArray.push( { data: data, type: ["draggableItem"], description: description  } );
		}

		// GENERATE DROP TARGET
		var dropTarget = new dojo.dnd.Source(
			this.dropTargetContainer,
			{
				//copyOnly: true,
				//selfAccept: false
				//,
				accept : [ "draggableItem" ]
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

			var project = {
				name : node.name,
				description : node.description
			};
			project.parentWidget = this;			

			var groupProjectRow = new plugins.admin.GroupProjectRow(project);
			node.innerHTML = '';
			node.appendChild(groupProjectRow.domNode);
		}


		// ADD SPECIAL DRAG AVATAR
		var groupProjectsObject = this;
		dropTarget.creator = function (item, hint)
		{

			// item CAN COME FROM THE DRAG SOURCE OR THE DROP TARGET ITSELF
			var name = item.projectname ? item.projectname : item.name;
			var description = item.projectdesc ? item.projectdesc : item.description;

			if ( name == null || name == '' )
			{
				return;
			}
			if ( description == null ) description = '';

			// ADD VALUES TO NODE SO THAT THEY GET PASSED TO this.addToGroup
			var node = dojo.doc.createElement("div");
			node.name = name;
			node.description = description;
			node.id = dojo.dnd.getUniqueId();
			node.className = "dojoDndItem";

			var project = {
				name : name,
				description : description
			};
			project.parentWidget = this;			

			var groupProjectRow = new plugins.admin.GroupProjectRow(project);
			node.innerHTML = '';
			node.appendChild(groupProjectRow.domNode);

			// ACCEPTABLE NODE DROPPED ONTO SELF --> ADD TO USER GROUP.
			if ( hint != 'avatar' )
			{
				groupProjectsObject.addToGroup(node.name, node.description, node.name);
			}

			return {node: node, data: item, type: ["draggableItem"]};
		};


		// ADD NODE IF DROPPED FROM OTHER SOURCE, DELETE IF DROPPED FROM SELF
		dojo.connect(dropTarget, "onDndDrop", function(source, nodes, copy, target)
		{
			// NODE DROPPED FROM SELF --> DELETE NODE
			if ( source == this && target != this )
			{
				for ( var i = 0; i < nodes.length; i++ )
				{
					var node = nodes[i];
					groupProjectsObject.removeFromGroup(node.name, node.description, node.location);
				}
			}			
			else
			{
				// DO NOTHING IF NODE WAS DROPPED FROM DRAG SOURCE
				// ( IT HAS ALREADY BEEN GENERATED BY dropTarget.creator() )
			}			

			// REMOVE DUPLICATE NODES
			var currentNodes = dropTarget.getAllNodes();
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
					// HACK TO AVOID THIS ERROR: node.parentNode is null
					try {
						node.parentNode.removeChild(node)

						// NB: THIS HAS NO EFFECT ON NODES IN dnd.Source
						//currentNodes.splice(i, 1); 
					}
					catch (e) {}
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
		groupObject.username = Agua.cookie('username');

		// ADD SOURCE OBJECT TO THE SOURCES IN THIS GROUP
		var groupName = this.groupCombo.getValue();

		var success = Agua.addProjectToGroup(groupName, groupObject);
		if ( success == false )
		{
			return;
		}

		// ADD THE SOURCE INTO THE groupprojects TABLE ON THE SERVER
		var data = new Object;
		data.name = name;
		data.description = description;
		data.location = location;
		data.groupname = groupName;
		data.type = "project";

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
		groupObject.username = Agua.cookie('username');

		if ( Agua.removeProjectFromGroup(groupName, groupObject) == false )
		{
			return;
		}

		// REMOVE THE SOURCE FROM THE groupprojects TABLE ON THE SERVER
		var data = new Object;
		data.name = String(name);
		data.description = description;
		data.location = location;
		data.groupname = groupName;
		data.type = "project";

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

	// setTrash
	//
	//	DELETE NODE IF DROPPED INTO TRASH. ACTUAL REMOVAL FROM THE
	//	DATA IS ACCOMPLISHED IN THE onDndDrop LISTENER OF THE SOURCE
	//
	setTrash : function ()
	{

		var trash = new dojo.dnd.Source(
			this.trashContainer,
			{
				accept : [ "draggableItem" ]
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

}); // plugins.admin.GroupProjects
