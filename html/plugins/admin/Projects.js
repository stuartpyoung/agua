dojo.provide("plugins.admin.Projects");

// ALLOW THE USER TO ADD, REMOVE AND MODIFY PROJECTS

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dijit.form.Button");
dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.Textarea");
dojo.require("dojo.parser");
dojo.require("dojo.dnd.Source");
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.admin.ProjectRow");

dojo.declare(
    "plugins.admin.Projects",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/projects.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	//addingProject STATE
	addingProject : false,

	// OR USE @import IN HTML TEMPLATE
	cssFiles : [ "plugins/admin/css/projects.css" ],

	// PARENT WIDGET
	parentWidget : null,

	constructor : function(args) {
		// GET INFO FROM ARGS
		this.parentWidget = args.parentWidget;
		this.projects = args.parentWidget.projects;

        // LOAD SORIA AND FILEPICKER CSS
        this.loadCSS();		
	},

	postCreate : function() {

		this.startup();
	},


	startup : function () {

		// COMPLETE CONSTRUCTION OF OBJECT
		this.inherited(arguments);	 

		// ADD TO TAB CONTAINER		
		this.tabContainer.addChild(this.projectsTab);
		this.tabContainer.selectChild(this.projectsTab);

		// SET DRAG PROJECT - LIST OF PROJECTS
		this.setDragSource();

		// SET NEW PROJECT FORM
		this.setNewProject();

		// SET TRASH DROP TARGET
		this.setTrash();
	},

	// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS
	reloadProjectTabs : function () {

		var tabPaneNames = ["plugins.admin.GroupProjects"];
		for ( var i in tabPaneNames )
		{
			if ( this.parentWidget.paneWidgets[tabPaneNames[i]] != null )
			{
				this.parentWidget.paneWidgets[tabPaneNames[i]].reload();
			}
		}
	},


	clearValue : function (textarea, value) {

		//if ( textarea.clicked == true ) return;
		//textarea.clicked = true;

		//if ( textarea.value.match(/^\s*/ + value + /\s*$/) )
		if ( textarea.value.match(value) )
		{
			textarea.value = '';
			//textarea.focus();
		}
	},

	setNewProject : function () {

		// FOCUS ON WIDGET TO FIX ITS WIDTH
		this.newName.focus();
		//this.addProjectButton.focus();

		// SET ADD PROJECT ONCLICK
		dojo.connect(this.addProjectButton, "onClick", dojo.hitch(this, "addProject"));

		// SET ONCLICK TO CANCEL DEFAULT TEXT
		dojo.connect(this.newName, "onclick", dojo.hitch(this, "clearValue", this.newName, 'Project'));
		dojo.connect(this.newName, "onfocus", dojo.hitch(this, "clearValue", this.newName, 'Project'));
		dojo.connect(this.newDescription, "onclick", dojo.hitch(this, "clearValue", this.newDescription, 'Project description'));
		dojo.connect(this.newDescription, "onfocus", dojo.hitch(this, "clearValue", this.newDescription, 'Project description'));
		dojo.connect(this.newNotes, "onclick", dojo.hitch(this, "clearValue", this.newNotes, 'Project notes'));
		dojo.connect(this.newNotes, "onfocus", dojo.hitch(this, "clearValue", this.newNotes, 'Project notes'));

		// SET NEW PROJECT LISTENER
		var projectsObject = this;
		dojo.connect(this.newName, "onkeypress", function(evt){
			var key = evt.charOrCode;
			if ( key == 13 )
			{				
				// SHIFT FOCUS TO NEXT INPUT
				projectsObject.newDescription.focus();
			}

			// REMOVE ANY LINE RETURNS
			this.value = this.value.replace(/\n/, '');
		});

		// SET NEW PROJECT LISTENER
		var projectsObject = this;
		dojo.connect(this.newDescription, "onkeypress", function(evt){
			var key = evt.charOrCode;
			if ( key == 13 )
			{
				// SHIFT FOCUS TO NEXT INPUT
				projectsObject.newNotes.focus();
			}

			// REMOVE ANY LINE RETURNS
			this.value = this.value.replace(/\n/, '');
		});

		// SET NEW PROJECT LISTENER
		var projectsObject = this;
		dojo.connect(this.newNotes, "onkeypress", function(evt){
			var key = evt.charOrCode;
			if ( key == 13 )
			{
				// SHIFT FOCUS TO SAVE BUTTON
				projectsObject.addProjectButton.focus();
			}

			// REMOVE ANY LINE RETURNS
			this.value = this.value.replace(/\n/, '');
		});

	},



	setDragSource : function () {

		// DELETE EXISTING TABLE CONTENT
		while ( this.dragSource.firstChild )
		{
			this.dragSource.removeChild(this.dragSource.firstChild);
		}

		var dataArray = new Array;
		var sourceArray = dojo.clone(Agua.projects);

		sourceArray = this.sortHasharray(sourceArray, 'name');

		// CHECK sourceArray IS NOT NULL OR EMPTY
		if ( sourceArray == null || sourceArray.length == 0 )
		{
			return;
		}

		// GENERATE dataArray TO INSERT INTO DND PROJECT TABLE
		for ( var j = 0; j < sourceArray.length; j++ )
		{
			var data = sourceArray[j];				
			data.toString = function () { return this.name; }
			dataArray.push( { data: data, type: ["draggableItem"] } );
		}

		// GENERATE DND PROJECT
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
			if ( node.description == null || node.description == '' )
			{
				node.description = "Project description";
				dojo.addClass(node, 'default');
			}
			node.notes = dataArray[k].data.notes;
			if ( node.notes == null || node.notes == '' )
			{
				node.notes = "Project notes";
				dojo.addClass(node, 'default');
			}

			var source = {
				name : node.name,
				description : node.description,
				notes : node.notes
			};

			source.parentWidget = this;

			var sourceRow = new plugins.admin.ProjectRow(source);

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


	editProjectRow : function (projectRow, node) {

		// RETURN IF ALREADY EDITING PROJECT ROW (I.E., MULTIPLE CLICKS)
		if ( this.editingProjectRow == true ) return;
		this.editingProjectRow = true;

		// REPLACE THE TD INNERHTML WITH A TEXTAREA
		var text = node.innerHTML;
		if ( text.match(/Project notes/) || text.match(/Project description/) )
			text = '';

		node.innerHTML = '';
		if ( text == null || ! text ) text = '';

		// RETURN IF THIS IS A DOUBLE-CLICK
		if ( text.match(/^<i/) ||
			text.match(/^<br/) ||
			text.match(/^<fieldset/) ||
			text.match(/^<textarea/) )
		{
			this.editingProjectRow = false;
			return;
		}

		// CREATE INPUT TEXT AREA
		var textarea = document.createElement('textarea');
		dojo.addClass(textarea, 'editProjectRow');
		node.appendChild(textarea);
		textarea.value = text;
		textarea.focus();

		// SET NEW PROJECT LISTENER
		var projectsObject = this;
		dojo.connect(textarea, "onkeypress", function(evt){

			// summary: handles keyboard events
			var key = evt.charOrCode;

			if ( key == 13 )
			{
				var newText = textarea.value;

				var project = new Object;
				project.name = projectRow.name.innerHTML;
				project.description = projectRow.description.innerHTML;
				project.notes = projectRow.notes.innerHTML;

				if ( project.description.match(/^<textarea/) )
					project.description = projectRow.description.firstChild.value;
				//if ( project.description.match(/^<i class/) )
				//	project.description =
				//	projectRow.description.firstChild.firstChild ?
				//	projectRow.description.firstChild.firstChild.value :
				//	projectRow.description.firstChild.innerHTML;
				if ( project.description == 'Project description' )
					project.description = '';

				if ( project.notes.match(/^<textarea/) )
					project.notes = projectRow.notes.firstChild.value;
				//if ( project.notes.match(/^<i class/) )
				//	project.notes =
				//	projectRow.notes.firstChild.firstChild ?
				//	projectRow.notes.firstChild.firstChild.value :
				//	projectRow.notes.firstChild.innerHTML;
				if ( project.notes == 'Project notes' )
					project.notes = '';



				// REMOVE WHITESPACE
				project.name = project.name.replace(/^[\s\n]+/, '');
				project.name = project.name.replace(/[\s\n]+$/, '');
				project.description = project.description.replace(/^[\s\n]+/g, '');
				project.description = project.description.replace(/[\s\n]+$/g, '');
				project.notes = project.notes.replace(/^[\s\n]+/, '');
				project.notes = project.notes.replace(/[\s\n]+$/, '');


				// REMOVE TEXTAREA
				node.removeChild(textarea);
				if ( dojo.hasClass(node, "description") && newText == '' )
					newText = "Project description";
				if ( dojo.hasClass(node, "notes") && newText == '' )
					newText = "Project notes";
				node.innerHTML = newText;

				// UNSET FLAG
				projectsObject.editingProjectRow = false;

				// SAVE NEW PROJECT TO REMOTE DATABASE
				projectsObject.saveProject(project);
			}
			else if (key == dojo.keys.ESCAPE)
			{
				projectsObject.editingProjectRow = false;

				// REMOVE TEXTAREA
				node.removeChild(textarea);
				if ( dojo.hasClass(node.parentNode, "description") && text == '' )
					text = "Project description";
				if ( dojo.hasClass(node.parentNode, "notes") && rext == '' )
					text = "Project notes";
				node.innerHTML = text;

			}
		});

	},



	deleteProject : function (name) {

		// CLEAN UP WHITESPACE
		name = name.replace(/\s+$/,'');
		name = name.replace(/^\s+/,'');

		var sourceObject = { name: name };

		// REMOVING PROJECT FROM Agua.projects
		var success = Agua.removeProject(sourceObject)

		// RESET THE PROJECTS TABLE
		this.setDragSource();

		var url = Agua.cgiUrl + "/agua?";

		// CREATE JSON QUERY
		var query = new Object;
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.mode = "deleteProject";
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

		// RELOAD RELEVANT DISPLAYS IN PROJECT-RELATED TABS
		this.reloadProjectTabs();

	}, // Projects.deleteProject


	addProject : function (source) {

		var name = this.newName.value;
		name = name.replace(/\s+/g, '');
		var description = this.newDescription.value;
		var notes = this.newNotes.value;

		// CHECK FOR VALID INPUTS
		if ( name == '' || name.match(/^\s*Project\s*$/) )
		{
			dojo.addClass(this.newName, 'invalid');
		}
		else{
			dojo.removeClass(this.newName, 'invalid');
		}
		if ( description == '' || description.match(/^\s*Description\s*$/) )
		{
			dojo.addClass(this.newDescription, 'invalid');
		}
		else{
			dojo.removeClass(this.newDescription, 'invalid');
		}
		if ( notes == '' || notes.match(/^\s*Notes\s*$/) )
		{
			dojo.addClass(this.newNotes, 'invalid');
		}
		else{
			dojo.removeClass(this.newNotes, 'invalid');
		}

		if ( name == '' || name.match(/^\s*Project\s*$/) 
			|| description == '' || description.match(/^\s*Description\s*$/)
			|| notes == '' || notes.match(/^\s*Notes\s*$/) )
		{
			return;
		}

		var sourceObject = { name: name, description: description, notes: notes };

		if ( Agua.isProject(sourceObject.name) )
		{
			dojo.removeClass(this.newName, 'invalid');
			dojo.addClass(this.newName, 'invalid');
			this.newName.focus();
			return;
		}

		// ADD PROJECT TO Agua.projects ARRAY
		Agua.addProject(sourceObject);

		// RESET PROJECT TABLE
		this.setDragSource();

		// SAVE PROJECT TO REMOTE DATABASE
		this._addProject(sourceObject);

		// RELOAD RELEVANT DISPLAYS IN PROJECT-RELATED TABS
		this.reloadProjectTabs();

	},	// Projects.addProject


	_addProject : function (source) {

		if ( this.addingProject == true )
		{
			return;
		}
		this.addingProject = true;

		// CLEAN UP WHITESPACE AND SUBSTITUTE NON-JSON SAFE CHARACTERS
		source.name = this.jsonSafe(source.name, 'toJson');
		source.description = this.jsonSafe(source.description, 'toJson');
		source.notes = this.jsonSafe(source.notes, 'toJson');

		var url = Agua.cgiUrl + "/agua?";

		// CREATE JSON QUERY
		var query = new Object;
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.mode = "addProject";
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

		this.addingProject = false;

	}, // Projects._addProject


	// UPDATE USER IN REMOTE DATABASE
	saveProject : function (savedObject) {

		if ( this.savingProject == true )
		{
			return;
		}
		this.savingProject = true;

		// CLEAN UP WHITESPACE AND SUBSTITUTE NON-JSON SAFE CHARACTERS
		savedObject.name = this.jsonSafe(savedObject.name, 'toJson');
		savedObject.description = this.jsonSafe(savedObject.description, 'toJson');
		savedObject.notes = this.jsonSafe(savedObject.notes, 'toJson');

		// REMOVE ORIGINAL project OBJECT FROM Agua.projects ARRAY
		// THEN ADD NEW PROJECT OBJECT TO Agua.projects ARRAY
		Agua.removeProject({ name: savedObject.name});
		Agua.addProject(savedObject);

		var url = Agua.cgiUrl + "agua";
		var query = savedObject;
		savedObject.username = Agua.cookie('username');
		savedObject.sessionId = Agua.cookie('sessionId');
		savedObject.mode = "saveProject";

		var returned = this.doPut({ url: url, query: query, sync: false });

		this.savingProject = false;

	}, // Projects.saveProject

	//	DELETE NODE IF DROPPED INTO TRASH. ACTUAL REMOVAL FROM THE
	//	DATA IS ACCOMPLISHED IN THE onDndDrop LISTENER OF THE PROJECT
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
						//node.parentNode.removeChild(node);

						thisObject.deleteProject(node.name);
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

}); // plugins.admin.Projects


