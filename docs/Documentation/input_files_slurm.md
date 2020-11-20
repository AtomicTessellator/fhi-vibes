# Slurm Submission

`vibes` can submit jobs on cluster using the `slurm` submission system. The configuration goes via the task input file by adding a `[slurm]` section. A template can be generated via

```
vibes template slurm >> your_input.in
```

??? info "Task input file with `[slurm]` section"

	```
	...
	[slurm]
	name:         test
	tag:          vibes
	mail_type:    all
	mail_address: your@mail.com
	nodes:        1
	cores:        32
	queue:        express
	timeout:      30
	```

## Submit a job

A job is submitted by using `vibes submit` instead of `vibes run`. For example, you would run `vibes submit md md.in` on a cluster to submit the calculation to the queue instead of running `vibes run md md.in` locally. `vibes submit` will submit the job to the cluster according to the specification it finds in the `[slurm]` section, see [below](#the-slurm-section).

The command will do the following:

- Write a `submit.sh` script according to the specifications found in the `[slurm]` section,
- submit the job to the queue via `sbatch submit.sh`, and
- log the time of submission and the job ID to a file called `.submit.log`.

## Restart a job

`vibes` supports autmatic restarts if the job will take longer than the available walltime, which is often the case during [_ab initio_ molecular dynamics simulations](../Tutorial/3_md_ab_initio.md). Restarts can be requested by adding

```
...
[restart]
command: the command to restart the calculation, e.g., `vibes submit md md.in`
```

to your task input file. This will run the command specified here shortly before the walltime is over.

## The `[slurm`] section

### `name`

`str`: The name of the job.

### `tag`

`str`: A common tag added to the job `name` to cluster sets of jobs.

### `mail_type`

`str`:  The [slurm mail type usually specified via `--mail-type`](https://slurm.schedmd.com/sbatch.html)

### `mail_address`

`str`: The mail address that slurm will send notifications to.

### `nodes`

`int`: The number of nodes to be used for the job (job dependent).

### `cores`

`int`: The number of cores per node to be used for the job (**machine dependent!**).

### `queue`

`str`: The name of the queue to be submitted to.

### `timeout`

`int`: Walltime for the job in minute. **Queue and cluster dependent!**