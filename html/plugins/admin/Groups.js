dojo.provide("plugins.admin.Groups");

// ALLOW THE USER TO ADD, REMOVE AND MODIFY GROUPS

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dijit.form.Button");
dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.Textarea");
dojo.require("dojo.parser");
dojo.require("dojo.dnd.Source");
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.admin.GroupRow");

dojo.declare(
    "plugins.admin.Groups",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/groups.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	//addingGroup STATE
	addingGroup : false,

	// OR USE @import IN HTML TEMPLATE
	cssFiles : [ "plugins/admin/css/groups.css" ],

	// PARENT WIDGET
	parentWidget : null,


	constructor : function(args)
	{
		// GET INFO FROM ARGS
		this.parentWidget = args.parentWidget;
		this.groups = args.parentWidget.groups;

        // LOAD CSS
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

		// ADD TO TAB CONTAINER		
		this.tabContainer.addChild(this.groupsTab);
		this.tabContainer.selectChild(this.groupsTab);

		// SET DRAG GROUP - LIST OF GROUPS
		this.setDragSource();

		// SET NEW GROUP FORM
		this.setForm();

		// SET TRASH DROP TARGET
		this.setTrash();

		// RESIZE ADMIN TAB ON window.resize
		dojo.connect(window,"resize",dojo.hitch(this, function(){
			this.groupsTab.resize();
			this.groupsTab.refresh();
		}));

	},


	// reloadGroupTabs
	//
	// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS
	//
	reloadGroupTabs : function ()
	{

		//for ( var tabPaneName in this.parentWidget.paneWidgets )
		//{
		//}

		var tabPaneNames = [ "plugins.admin.GroupSources", "plugins.admin.GroupUsers", "plugins.admin.GroupProjects" ];
		//var tabPaneNames = [ "plugins.admin.GroupSources", "plugins.admin.GroupUsers" ];
		for ( var i in tabPaneNames )
		{
			if ( this.parentWidget.paneWidgets[tabPaneNames[i]] != null )
			{
				this.parentWidget.paneWidgets[tabPaneNames[i]].reload();
			}
		}
	},

	// setForm
	//
	// SET 'ADD NEW GROUP' FORM
	//
	setForm : function ()
	{

		// FOCUS ON WIDGET TO FIX ITS WIDTH
		this.newName.focus();
		//this.addGroupButton.focus();

		// SET ADD GROUP ONCLICK
		dojo.connect(this.addGroupButton, "onClick", dojo.hitch(this, "addGroup"));

		// SET ONCLICK TO CANCEL DEFAULT TEXT
		dojo.connect(this.newName, "onclick", dojo.hitch(this, "clearValue", this.newName, 'Group'));
		dojo.connect(this.newName, "onfocus", dojo.hitch(this, "clearValue", this.newName, 'Group'));
		dojo.connect(this.newDescription, "onclick", dojo.hitch(this, "clearValue", this.newDescription, 'Description'));
		dojo.connect(this.newDescription, "onfocus", dojo.hitch(this, "clearValue", this.newDescription, 'Description'));
		dojo.connect(this.newNotes, "onclick", dojo.hitch(this, "clearValue", this.newNotes, 'Group notes'));
		dojo.connect(this.newNotes, "onfocus", dojo.hitch(this, "clearValue", this.newNotes, 'Group notes'));

		// SET NEW PROJECT LISTENER
		var groupsObject = this;
		dojo.connect(this.newName, "onkeypress", function(evt){
			var key = evt.charOrCode;
			if ( key == 13 )
			{				
				// SHIFT FOCUS TO NEXT INPUT
				groupsObject.newDescription.focus();
			}

			// REMOVE ANY LINE RETURNS
			this.value = this.value.replace(/\n/, '');
		});

		// SET NEW PROJECT LISTENER
		var groupsObject = this;
		dojo.connect(this.newDescription, "onkeypress", function(evt){
			var key = evt.charOrCode;
			if ( key == 13 )
			{
				// SHIFT FOCUS TO NEXT INPUT
				groupsObject.newNotes.focus();
			}

			// REMOVE ANY LINE RETURNS
			this.value = this.value.replace(/\n/, '');
		});

		// SET NEW PROJECT LISTENER
		var groupsObject = this;
		dojo.connect(this.newNotes, "onkeypress", function(evt){
			var key = evt.charOrCode;
			if ( key == 13 )
			{
				// SHIFT FOCUS TO SAVE BUTTON
				groupsObject.addGroupButton.focus();
			}

			// REMOVE ANY LINE RETURNS
			this.value = this.value.replace(/\n/, '');
		});

	},



	setDragSource : function ()
	{

		// DELETE EXISTING TABLE CONTENT
		while ( this.dragSource.firstChild )
		{
			this.dragSource.removeChild(this.dragSource.firstChild);
		}

		var dataArray = new Array;
		var sourceArray = dojo.clone(Agua.groups);

		sourceArray = this.sortHasharray(sourceArray, 'name');

		// CHECK sourceArray IS NOT NULL OR EMPTY
		if ( sourceArray == null || sourceArray.length == 0 )
		{
			return;
		}

		// GENERATE dataArray TO INSERT INTO DND GROUP TABLE
		for ( var j = 0; j < sourceArray.length; j++ )
		{
			var data = sourceArray[j];				
			data.toString = function () { return this.name; }
			dataArray.push( { data: data, type: ["draggableItem"] } );
		}

		// GENERATE DND GROUP
		var dragSource = new dojo.dnd.Source(
			this.dragSource,
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


			// SET NODE name AND description
			node.name = dataArray[k].data.name;
			node.description = dataArray[k].data.description;
			node.notes = dataArray[k].data.notes;
			if ( node.notes == null || node.notes == '' )
			{
				node.notes = "<i class='default'>Group notes</i>";
			}

			var source = {
				name : node.name,
				description : node.description,
				notes : node.notes
			};

			source.parentWidget = this;

			var sourceRow = new plugins.admin.GroupRow(source);

			//node.innerHTML = sourceRow.domNode.innerHTML;
			node.innerHTML = '';
			node.appendChild(sourceRow.domNode);
			//document.body.appendChild(sourceRow.domNode);


		}

		var sourceObject = this;
		dragSource.creator = function (item, hint)
		{

			var node = dojo.doc.createElement("div");
			node.name = item.name;
			node.description = item.description;
			node.notes = item.notes;
			node.id = dojo.dnd.getUniqueId();
			node.className = "dojoDndItem";


			// SET FANCY FORMAT IN NODE INNERHTML
			node.innerHTML = "<table> <tr><td><strong style='color: darkred'>" + item.name + "</strong></td></tr><tr><td> " + item.description + "</td></tr></table>";

			return {node: node, data: item, type: ["text"]};
		};
	},


	editGroupRow : function (groupRow, node)
	{

		// RETURN IF ALREADY EDITING GROUP ROW (I.E., MULTIPLE CLICKS)
		if ( this.editingGroupRow == true )
		{
			return;
		}
		this.editingGroupRow = true;

		// REPLACE THE TD INNERHTML WITH A TEXTAREA
		var text = node.innerHTML;
		if ( text.match(/Group notes/) )
		{
			text = '';
		}

		node.innerHTML = '';
		if ( text == null || ! text ) text = '';

		// RETURN IF THIS IS A DOUBLE-CLICK
		if ( text.match(/^<i/) ||
			text.match(/^<br/) ||
			text.match(/^<fieldset/) ||
			text.match(/^<textarea/) )
		{
			this.editingGroupRow = false;
			return;
		}


		// CREATE INPUT TEXT AREA
		var textarea = document.createElement('textarea');
		dojo.addClass(textarea, 'editGroupRow');
		node.appendChild(textarea);
		textarea.value = text;
		textarea.focus();

		// SET NEW PROJECT LISTENER
		var groupsObject = this;
		dojo.connect(textarea, "onkeypress", function(evt){

			// summary: handles keyboard events
			var key = evt.charOrCode;

			if ( key == 13 )
			{
				var newText = textarea.value;

				// REMOVE TEXTAREA
				if ( node.firstChild == textarea )	node.removeChild(textarea);

				var newGroup = new Object;
				newGroup.name = groupRow.name.innerHTML;
				newGroup.description = groupRow.description.innerHTML;
				newGroup.notes = groupRow.notes.innerHTML;

				if ( newGroup.description.match(/^<textarea/) )
					newGroup.description = groupRow.description.firstChild.value;

				if ( newGroup.notes.match(/^<textarea/) )
					newGroup.notes = groupRow.notes.firstChild.value;

				// REMOVE WHITESPACE
				newGroup.name = newGroup.name.replace(/^\s+/, '');
				newGroup.name = newGroup.name.replace(/\s+$/, '');
				newGroup.description = newGroup.description.replace(/^\s+/, '');
				newGroup.description = newGroup.description.replace(/\s+$/, '');
				newGroup.notes = newGroup.notes.replace(/^\s+/, '');
				newGroup.notes = newGroup.notes.replace(/\s+$/, '');



				if ( newGroup.name != '' )
				{
					// REMOVE ORIGINAL newGroup OBJECT FROM Agua.groups ARRAY
					// THEN ADD NEW GROUP OBJECT TO Agua.groups ARRAY
					Agua.removeGroup({ name: newGroup.name});
					Agua.addGroup(newGroup);

					// REMOVE TEXTAREA
					node.removeChild(textarea);
					node.innerHTML = newText;

					// REDO GROUP TABLE
					groupsObject.setDragSource();

					// SAVE NEW GROUP TO REMOTE DATABASE
					groupsObject.saveGroup(newGroup);

					groupsObject.editingGroupRow = false;

				}
				else
				{
					groupsObject.editingGroupRow = false;
				}
			}
			else if (key == dojo.keys.ESCAPE)
			{
				groupsObject.editingGroupRow = false;

				// REMOVE TEXTAREA
				node.removeChild(textarea);
				if ( text == null || text == '' )
				{
					text = "<i class='default'>Group notes</i>";
				}
				node.innerHTML = text;

			}
		});

	},



	deleteGroup : function (name)
	{

		// CLEAN UP WHITESPACE
		name = name.replace(/\s+$/,'');
		name = name.replace(/^\s+/,'');

		var sourceObject = { name: name };

		// REMOVING GROUP FROM Agua.groups
		Agua.removeGroup(sourceObject)

		// RESET THE GROUPS TABLE
		this.setDragSource();

		var url = Agua.cgiUrl + "/agua?";

		// CREATE JSON QUERY
		var query = new Object;
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.mode = "deleteGroup";
		query.data = sourceObject;

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

		// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS
		this.reloadGroupTabs();

	}, // Groups.deleteGroup


	addGroup : function (source)
	{

		var name = this.newName.value;
		name = name.replace(/\s+/g, '');
		var description = this.newDescription.value;
		var notes = this.newNotes.value;

		// CHECK FOR VALID INPUTS
		/* name */
		if ( name == '' || name.match(/^\s*Group\s*$/) )
		{
			dojo.addClass(this.newName, 'invalid');
		}
		else{
			dojo.removeClass(this.newName, 'invalid');
		}

		/* description */
		//if ( description == '' || description.match(/^\s*Description\s*$/) )
		//{
		//	dojo.addClass(this.newDescription, 'invalid');
		//}
		//else{
		//	dojo.removeClass(this.newDescription, 'invalid');
		//}

		/* notes */
		//if ( notes == '' || notes.match(/^\s*Notes\s*$/) )
		//{
		//	dojo.addClass(this.newNotes, 'invalid');
		//}
		//else{
		//	dojo.removeClass(this.newNotes, 'invalid');
		//}

		if ( name == '' || name.match(/^\s*Group\s*$/) 
			|| description == '' || description.match(/^\s*Description\s*$/) )
			//|| notes == '' || notes.match(/^\s*Notes\s*$/) )
		{
			return;
		}

		var sourceObject = { name: name, description: description, notes: notes };

		if ( Agua.isGroup(sourceObject.name) )
		{
			dojo.removeClass(this.newName, 'invalid');
			dojo.addClass(this.newName, 'invalid');
			this.newName.focus();
			return;
		}

		// ADD GROUP TO Agua.groups ARRAY
		Agua.addGroup(sourceObject);

		// RESET GROUP TABLE
		this.setDragSource();

		// SAVE GROUP TO REMOTE DATABASE
		this.saveGroup(sourceObject);

		// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS
		this.reloadGroupTabs();

	},	// Groups.addGroup



	saveGroup : function (source)
	{

		if ( this.savingGroup == true )
		{
			return;
		}
		this.savingGroup = true;

		// CLEAN UP WHITESPACE AND SUBSTITUTE NON-JSON SAFE CHARACTERS
		source.name = this.jsonSafe(source.name, 'toJson');
		source.description = this.jsonSafe(source.description, 'toJson');
		source.notes = this.jsonSafe(source.notes, 'toJson');

		var url = Agua.cgiUrl + "/agua?";

		// CREATE JSON QUERY
		var query = new Object;
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.mode = "saveGroup";
		query.data = source;

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

		this.savingGroup = false;

	}, // Groups.saveGroup




	// setTrash
	//
	//	DELETE NODE IF DROPPED INTO TRASH. ACTUAL REMOVAL FROM THE
	//	DATA IS ACCOMPLISHED IN THE onDndDrop LISTENER OF THE GROUP
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
		var thisObject = this;
		dojo.connect(trash, "onDndDrop", function(source, nodes, copy, target){
			// NODE DROPPED ON SELF --> DELETE THE NODE
			if ( target == this )
			{

				for ( var i = 0; i < nodes.length; i++ )
				{
					var node = nodes[i];

					// HACK TO AVOID THIS ERROR: node.parentNode is null
					try {
						node.parentNode.removeChild(node);

						thisObject.deleteGroup(node.name);
					}
					catch (e) {
					}
				}

				// DELETE EXISTING TABLE CONTENT
				while ( thisObject.trashContainer.childNodes.length > 2 )
				{

					thisObject.trashContainer.removeChild(thisObject.trashContainer.childNodes[2]);
				}

			}

		});
	}

}); // plugins.admin.Groups

