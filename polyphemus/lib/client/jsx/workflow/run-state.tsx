export const RUN_NEVER = -1;
export const RUN_INTERVAL = 1;

export const getRunState = (run:number) => (run > 0 ? RUN_INTERVAL : run);
export const getRunIntervalTime = (run:number) => (run > 0 ? run : 0);

export const formatTime = (time:Date|string) => Intl.DateTimeFormat('en-US', {
   year: 'numeric', month: 'short', day: 'numeric',
   hour: 'numeric', minute: 'numeric'
}).format(new Date(time));

export const runDesc = (run:number) => {
  if (run == RUN_NEVER) return 'never';
  return `every ${run}s`;
}

export const runTime = (ran_at:string, run:number) => {
  if (run == RUN_NEVER) return 'never';
  if (!ran_at) return 'pending';

  const ranAt = new Date(ran_at);
  const runAt = new Date( ranAt.getTime() + run * 1000 );
  const now = new Date();

  if (runAt < now) return 'pending';

  return formatTime( runAt );
};

