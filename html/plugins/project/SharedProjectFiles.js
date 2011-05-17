
dojo.provide("plugins.project.SharedProjectFiles");

// **********************
// *** LOAD SHARED PROJECTS  ****
// **********************

// LATER FOR MENU: HOW TO DYNAMICALLY
// ENABLE / DISABLE MENU ITEM
// attr('disabled', bool) 

// DISPLAY THE USER'S PROJECTS DIRECTORY AND ALLOW
// THE USER TO BROWSE FILES AND MANIPULATE WORKFLOW
// FOLDERS AND FILES

// INHERITS
dojo.require("plugins.project.ProjectFiles");

dojo.declare( "plugins.project.SharedProjectFiles",
	[ plugins.project.ProjectFiles ],
{

	// *** LOAD FILESYSTEM ****
	load : function ()
	{


		var usernames = Agua.getSharedUsernames();
		for ( var j = 0; j < usernames.length; j++ )
		{
			var shares = Agua.getSharedProjectsByUsername(usernames[j]);

			for ( var i = 0; i < shares.length; i++ )
			{

				var name = shares[i].project;
				var description = shares[i].description;
				var owner = shares[i].owner;
				var groupname = shares[i].groupname;
				if ( ! description ) { description = '' };

				var titlePane = this.createTitlePane(
				{
					owner: owner,
					name: name,
					description: description,
					open: this.open
				});


				var thisStore = this.createStore(shares[i]);


				// GENERATE NEW FileDrag OBJECT
				var fileDrag = new plugins.files.FileDrag(
					{
						style: "height: auto; width: 100%; minHeight: 50px;",
						store: thisStore,
						fileMenu: this.fileMenu,
						folderMenu: this.folderMenu,
						workflowMenu: this.workflowMenu,
						owner: owner,
						parentWidget: this
					}
				);

				// SET PATH FOR THIS SHARE
				fileDrag.path = name;                    

				// START UP FileDrag
				fileDrag.startup();

				// ADD shareDrag TO TITLE PANE
				titlePane.containerNode.appendChild(fileDrag.domNode);

	//break;

			} // shares

		} // usernames
	},


	// CREATE URL FOR STORE
	createUrl : function (directory)
	{

		// SET URL
		var url = Agua.cgiUrl;
		url += "project.cgi?";
		url += "mode=fileSystem";
		url += "&sessionId=" + Agua.cookie('sessionId');
		url += "&requestor=" + Agua.cookie('username');
		url += "&username=" + directory.owner;
		url += "&groupname=" + directory.groupname;

		return url;		
	},

	// CREATE STORE FOR FILE DRAG
	createStore : function (directory)
	{

		// SET URL
		var url = this.createUrl(directory);

		// CREATE STORE
		var thisStore = new dojox.data.FileStore(
			{
				//id: paneNodeId + "-fileStore",
				url: url,
				pathAsQueryParam: true
			}
		);

		// SET FILE STORE path TO project
		thisStore.preamble = function()
		{
			this.store.path = this.arguments[0].path;                        
		};

		return thisStore;		
	}
});
