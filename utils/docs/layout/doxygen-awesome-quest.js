
/* injects the QuEST version number below the logo in the 
 * left nav-bar; this requires this script is included in
 * header.html, and the version number style is customised
 * in *-quest.css. This snippet was taken from Zephyr's
 * github (github.com/zephyrproject-rtos/zephyr) at:
 * zephyr/doc/_doxygen/doxygen-awesome-zephyr.js */

window.addEventListener('DOMContentLoaded', (event) => {
    let version = document.getElementById('projectnumber').innerText
    let titleTable = document.querySelector('#titlearea table');
    let cell = titleTable.insertRow(1).insertCell(0);
    cell.innerHTML = '<div id="projectversion">' + version + '</div>';
});
