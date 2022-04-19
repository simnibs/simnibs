.. _gallery:

Gallery
=========================
.. raw:: html

  <embed>
	<style>
	  .gallerytable {
		position: relative;
		width: 100%;
		}
	  .gallerytable > table {
		width: 100%;
		table-layout: fixed;
		}
	  .gallerytable > table > tbody > tr > td {
		background-size: contain;
		background-repeat: no-repeat;
		background-position: center center;
		text-align: center;
		vertical-align: bottom;
		overflow: hidden;
		text-overflow: ellipsis;
	   }
	  .gallerytable > table > tbody > tr > td > a {
		display: block;
		height: 100%;
		width: 100%;
		text-decoration: none;
	}
	</style>
	<div class="gallerytable"; id="gallerytable";></div>
	<script type="text/javascript" src=_static/gallery/list.js></script>
	<script>
		window.onload = function() {
			Nentries=filelist.length
			Ncolumn=2
			Nrow=Math.ceil(Nentries/Ncolumn)
			let table = document.createElement('table');
			entry=0
			for (let row = 0; row < Nrow; row++) {
			  table.insertRow();
			  table.insertRow();
			  for (let column = 0; column < Ncolumn; column++) {
				if (entry > Nentries - 1) { break; }
				let newCell = table.rows[table.rows.length - 2].insertCell();
				newCell.textContent = description[entry];
				newCell.style="padding-top:2vh;"
				let newCell2 = table.rows[table.rows.length - 1].insertCell();				
				newCell2.style="background-image:url(_static/gallery/"+filelist[entry]+");height:50vh;"
				if ( links[entry] ) {
					let aTag = document.createElement('a');
					aTag.setAttribute('href',links[entry]);
					newCell2.appendChild(aTag);
					}
				entry++
			  }
			}
			document.getElementById("gallerytable").appendChild(table);
		}
	</script>
  </embed>

